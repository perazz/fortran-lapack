!> Hessenberg reduction: balancing, back-transformation, orthogonal factor generation
module la_lapack_eigv_gen_hess
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_tri
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_scalar
     use la_lapack_householder_reflectors
     use la_lapack_orthogonal_factors_qr
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sgebak
     public :: la_sorghr
     public :: la_sormhr
     public :: la_sgebal
     public :: la_sgehd2
     public :: la_slahr2
     public :: la_sgehrd
     public :: la_dgebak
     public :: la_dorghr
     public :: la_dormhr
     public :: la_dgebal
     public :: la_dgehd2
     public :: la_dlahr2
     public :: la_dgehrd
#ifdef LA_WITH_XDP
     public :: la_xgebak
     public :: la_xorghr
     public :: la_xormhr
     public :: la_xgebal
     public :: la_xgehd2
     public :: la_xlahr2
     public :: la_xgehrd
#endif
#ifdef LA_WITH_QP
     public :: la_qgebak
     public :: la_qorghr
     public :: la_qormhr
     public :: la_qgebal
     public :: la_qgehd2
     public :: la_qlahr2
     public :: la_qgehrd
#endif
     public :: la_cgebak
     public :: la_cgebal
     public :: la_cgehd2
     public :: la_clahr2
     public :: la_cunghr
     public :: la_cunmhr
     public :: la_cgehrd
     public :: la_zgebak
     public :: la_zgebal
     public :: la_zgehd2
     public :: la_zlahr2
     public :: la_zunghr
     public :: la_zunmhr
     public :: la_zgehrd
#ifdef LA_WITH_XDP
     public :: la_ygebak
     public :: la_ygebal
     public :: la_ygehd2
     public :: la_ylahr2
     public :: la_yunghr
     public :: la_yunmhr
     public :: la_ygehrd
#endif
#ifdef LA_WITH_QP
     public :: la_wgebak
     public :: la_wgebal
     public :: la_wgehd2
     public :: la_wlahr2
     public :: la_wunghr
     public :: la_wunmhr
     public :: la_wgehrd
#endif

     contains

     !> SGEBAK: forms the right or left eigenvectors of a real general matrix
     !> by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by SGEBAL.

     pure subroutine la_sgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(inout) :: v(ldv,*)
           real(sp),intent(in) :: scale(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(sp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_sscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_sscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_sswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_sswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_sgebak
     !> DGEBAK: forms the right or left eigenvectors of a real general matrix
     !> by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by DGEBAL.

     pure subroutine la_dgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(in) :: scale(*)
           real(dp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(dp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_dscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_dscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_dswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_dswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_dgebak
#ifdef LA_WITH_XDP
     !> XGEBAK: forms the right or left eigenvectors of a real general matrix
     !> by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by XGEBAL.

     pure subroutine la_xgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(in) :: scale(*)
           real(xdp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(xdp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_xscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_xscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_xswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_xswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_xgebak
#endif
#ifdef LA_WITH_QP
     !> QGEBAK: forms the right or left eigenvectors of a real general matrix
     !> by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by QGEBAL.

     pure subroutine la_qgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(in) :: scale(*)
           real(qp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(qp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_qscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_qscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_qswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_qswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_qgebak
#endif

     !> SORGHR: generates a real orthogonal matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> SGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_sorghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'SORGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SORGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors one
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = zero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = zero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_sorgqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_sorghr
     !> DORGHR: generates a real orthogonal matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> DGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_dorghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'DORGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DORGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors one
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = zero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = zero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_dorgqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_dorghr
#ifdef LA_WITH_XDP
     !> XORGHR: generates a real orthogonal matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> XGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_xorghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'XORGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XORGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors one
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = zero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = zero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_xorgqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_xorghr
#endif
#ifdef LA_WITH_QP
     !> QORGHR: generates a real orthogonal matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> QGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_qorghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'QORGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QORGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors one
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = zero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = zero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = zero
              end do
              a(j,j) = one
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_qorgqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_qorghr
#endif

     !> SORMHR: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by SGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_sormhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),c(ldc,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'SORMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'SORMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SORMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_sormqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_sormhr
     !> DORMHR: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by DGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_dormhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),c(ldc,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'DORMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'DORMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DORMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_dormqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_dormhr
#ifdef LA_WITH_XDP
     !> XORMHR: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by XGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_xormhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'XORMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'XORMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XORMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_xormqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_xormhr
#endif
#ifdef LA_WITH_QP
     !> QORMHR: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by QGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_qormhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),c(ldc,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'QORMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'QORMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QORMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_qormqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_qormhr
#endif

     !> SGEBAL: balances a general real matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_sgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: scale(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sclfac = 2.0e+0_sp
           real(sp),parameter :: factor = 0.95e+0_sp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(sp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_sswap(l,a(1,j),1,a(1,m),1)
           call la_sswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (a(j,i) /= zero) cycle loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (a(i,j) /= zero) cycle loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_slamch('S')/la_slamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_snrm2(l - k + 1,a(k,i),1)
              r = la_snrm2(l - k + 1,a(i,k),lda)
              ica = la_isamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_isamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
                 if (la_sisnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('SGEBAL',-info)
                 return
              end if
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_sscal(n - k + 1,g,a(i,k),lda)
              call la_sscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_sgebal
     !> DGEBAL: balances a general real matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_dgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: scale(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sclfac = 2.0e+0_dp
           real(dp),parameter :: factor = 0.95e+0_dp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(dp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_dswap(l,a(1,j),1,a(1,m),1)
           call la_dswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (a(j,i) /= zero) cycle loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (a(i,j) /= zero) cycle loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_dlamch('S')/la_dlamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_dnrm2(l - k + 1,a(k,i),1)
              r = la_dnrm2(l - k + 1,a(i,k),lda)
              ica = la_idamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_idamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_disnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('DGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_dscal(n - k + 1,g,a(i,k),lda)
              call la_dscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_dgebal
#ifdef LA_WITH_XDP
     !> XGEBAL: balances a general real matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_xgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: scale(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sclfac = 2.0e+0_xdp
           real(xdp),parameter :: factor = 0.95e+0_xdp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(xdp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('XGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_xswap(l,a(1,j),1,a(1,m),1)
           call la_xswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (a(j,i) /= zero) cycle loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (a(i,j) /= zero) cycle loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_xlamch('S')/la_xlamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_xnrm2(l - k + 1,a(k,i),1)
              r = la_xnrm2(l - k + 1,a(i,k),lda)
              ica = la_ixamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_ixamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_xisnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('XGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_xscal(n - k + 1,g,a(i,k),lda)
              call la_xscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_xgebal
#endif
#ifdef LA_WITH_QP
     !> QGEBAL: balances a general real matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_qgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: scale(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sclfac = 2.0e+0_qp
           real(qp),parameter :: factor = 0.95e+0_qp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(qp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_qswap(l,a(1,j),1,a(1,m),1)
           call la_qswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (a(j,i) /= zero) cycle loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (a(i,j) /= zero) cycle loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_qlamch('S')/la_qlamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_qnrm2(l - k + 1,a(k,i),1)
              r = la_qnrm2(l - k + 1,a(i,k),lda)
              ica = la_iqamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_iqamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_qisnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('QGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_qscal(n - k + 1,g,a(i,k),lda)
              call la_qscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_qgebal
#endif

     !> SGEHD2: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_sgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('SGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              call la_slarfg(ihi - i,a(i + 1,i),a(min(i + 2,n),i),1,tau(i))
              aii = a(i + 1,i)
              a(i + 1,i) = one
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_slarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i) to a(i+1:ihi,i+1:n) from the left
              call la_slarf('LEFT',ihi - i,n - i,a(i + 1,i),1,tau(i),a(i + 1,i + 1),lda, &
                        work)
              a(i + 1,i) = aii
           end do
           return
     end subroutine la_sgehd2
     !> DGEHD2: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_dgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              call la_dlarfg(ihi - i,a(i + 1,i),a(min(i + 2,n),i),1,tau(i))
              aii = a(i + 1,i)
              a(i + 1,i) = one
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_dlarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i) to a(i+1:ihi,i+1:n) from the left
              call la_dlarf('LEFT',ihi - i,n - i,a(i + 1,i),1,tau(i),a(i + 1,i + 1),lda, &
                        work)
              a(i + 1,i) = aii
           end do
           return
     end subroutine la_dgehd2
#ifdef LA_WITH_XDP
     !> XGEHD2: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_xgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(xdp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('XGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              call la_xlarfg(ihi - i,a(i + 1,i),a(min(i + 2,n),i),1,tau(i))
              aii = a(i + 1,i)
              a(i + 1,i) = one
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_xlarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i) to a(i+1:ihi,i+1:n) from the left
              call la_xlarf('LEFT',ihi - i,n - i,a(i + 1,i),1,tau(i),a(i + 1,i + 1),lda, &
                        work)
              a(i + 1,i) = aii
           end do
           return
     end subroutine la_xgehd2
#endif
#ifdef LA_WITH_QP
     !> QGEHD2: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_qgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              call la_qlarfg(ihi - i,a(i + 1,i),a(min(i + 2,n),i),1,tau(i))
              aii = a(i + 1,i)
              a(i + 1,i) = one
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_qlarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i) to a(i+1:ihi,i+1:n) from the left
              call la_qlarf('LEFT',ihi - i,n - i,a(i + 1,i),1,tau(i),a(i + 1,i + 1),lda, &
                        work)
              a(i + 1,i) = aii
           end do
           return
     end subroutine la_qgehd2
#endif

     !> SLAHR2: reduces the first NB columns of A real general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an orthogonal similarity transformation
     !> Q**T * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by SGEHRD.

     pure subroutine la_slahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**t
                 call la_sgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,one,a(k + 1,i),1)
                 ! apply i - v * t**t * v**t to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**t * b1
                 call la_scopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_strmv('LOWER','TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1,nb), &
                            1)
                 ! w := w + v2**t * b2
                 call la_sgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i), &
                           1,one,t(1,nb),1)
                 ! w := t**t * w
                 call la_strmv('UPPER','TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,nb),1)

                 ! b2 := b2 - v2*w
                 call la_sgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-one,a(k + i,1),lda,t(1,nb) &
                           ,1,one,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_strmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_saxpy(i - 1,-one,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_slarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = one
              ! compute  y(k+1:n,i)
              call la_sgemv('NO TRANSPOSE',n - k,n - k - i + 1,one,a(k + 1,i + 1),lda,a(k + i,i), &
                         1,zero,y(k + 1,i),1)
              call la_sgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i),1, &
                        zero,t(1,i),1)
              call la_sgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,t(1,i),1, &
                        one,y(k + 1,i),1)
              call la_sscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_sscal(i - 1,-tau(i),t(1,i),1)
              call la_strmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_slacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_strmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,one,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,one,a(1, &
                     2 + nb),lda,a(k + 1 + nb,1),lda,one,y,ldy)
           call la_strmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,one,t,ldt,y, &
                     ldy)
           return
     end subroutine la_slahr2
     !> DLAHR2: reduces the first NB columns of A real general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an orthogonal similarity transformation
     !> Q**T * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by DGEHRD.

     pure subroutine la_dlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**t
                 call la_dgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,one,a(k + 1,i),1)
                 ! apply i - v * t**t * v**t to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**t * b1
                 call la_dcopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_dtrmv('LOWER','TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1,nb), &
                            1)
                 ! w := w + v2**t * b2
                 call la_dgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i), &
                           1,one,t(1,nb),1)
                 ! w := t**t * w
                 call la_dtrmv('UPPER','TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,nb),1)

                 ! b2 := b2 - v2*w
                 call la_dgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-one,a(k + i,1),lda,t(1,nb) &
                           ,1,one,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_dtrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_daxpy(i - 1,-one,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_dlarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = one
              ! compute  y(k+1:n,i)
              call la_dgemv('NO TRANSPOSE',n - k,n - k - i + 1,one,a(k + 1,i + 1),lda,a(k + i,i), &
                         1,zero,y(k + 1,i),1)
              call la_dgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i),1, &
                        zero,t(1,i),1)
              call la_dgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,t(1,i),1, &
                        one,y(k + 1,i),1)
              call la_dscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_dscal(i - 1,-tau(i),t(1,i),1)
              call la_dtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_dlacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_dtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,one,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,one,a(1, &
                     2 + nb),lda,a(k + 1 + nb,1),lda,one,y,ldy)
           call la_dtrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,one,t,ldt,y, &
                     ldy)
           return
     end subroutine la_dlahr2
#ifdef LA_WITH_XDP
     !> XLAHR2: reduces the first NB columns of A real general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an orthogonal similarity transformation
     !> Q**T * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by XGEHRD.

     pure subroutine la_xlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(xdp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**t
                 call la_xgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,one,a(k + 1,i),1)
                 ! apply i - v * t**t * v**t to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**t * b1
                 call la_xcopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_xtrmv('LOWER','TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1,nb), &
                            1)
                 ! w := w + v2**t * b2
                 call la_xgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i), &
                           1,one,t(1,nb),1)
                 ! w := t**t * w
                 call la_xtrmv('UPPER','TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,nb),1)

                 ! b2 := b2 - v2*w
                 call la_xgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-one,a(k + i,1),lda,t(1,nb) &
                           ,1,one,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_xtrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_xaxpy(i - 1,-one,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_xlarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = one
              ! compute  y(k+1:n,i)
              call la_xgemv('NO TRANSPOSE',n - k,n - k - i + 1,one,a(k + 1,i + 1),lda,a(k + i,i), &
                         1,zero,y(k + 1,i),1)
              call la_xgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i),1, &
                        zero,t(1,i),1)
              call la_xgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,t(1,i),1, &
                        one,y(k + 1,i),1)
              call la_xscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_xscal(i - 1,-tau(i),t(1,i),1)
              call la_xtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_xlacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_xtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,one,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,one,a(1, &
                     2 + nb),lda,a(k + 1 + nb,1),lda,one,y,ldy)
           call la_xtrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,one,t,ldt,y, &
                     ldy)
           return
     end subroutine la_xlahr2
#endif
#ifdef LA_WITH_QP
     !> QLAHR2: reduces the first NB columns of A real general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an orthogonal similarity transformation
     !> Q**T * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by QGEHRD.

     pure subroutine la_qlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**t
                 call la_qgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,one,a(k + 1,i),1)
                 ! apply i - v * t**t * v**t to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**t * b1
                 call la_qcopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_qtrmv('LOWER','TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1,nb), &
                            1)
                 ! w := w + v2**t * b2
                 call la_qgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i), &
                           1,one,t(1,nb),1)
                 ! w := t**t * w
                 call la_qtrmv('UPPER','TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,nb),1)

                 ! b2 := b2 - v2*w
                 call la_qgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-one,a(k + i,1),lda,t(1,nb) &
                           ,1,one,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_qtrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_qaxpy(i - 1,-one,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_qlarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = one
              ! compute  y(k+1:n,i)
              call la_qgemv('NO TRANSPOSE',n - k,n - k - i + 1,one,a(k + 1,i + 1),lda,a(k + i,i), &
                         1,zero,y(k + 1,i),1)
              call la_qgemv('TRANSPOSE',n - k - i + 1,i - 1,one,a(k + i,1),lda,a(k + i,i),1, &
                        zero,t(1,i),1)
              call la_qgemv('NO TRANSPOSE',n - k,i - 1,-one,y(k + 1,1),ldy,t(1,i),1, &
                        one,y(k + 1,i),1)
              call la_qscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_qscal(i - 1,-tau(i),t(1,i),1)
              call la_qtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_qlacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_qtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,one,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,one,a(1, &
                     2 + nb),lda,a(k + 1 + nb,1),lda,one,y,ldy)
           call la_qtrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,one,t,ldt,y, &
                     ldy)
           return
     end subroutine la_qlahr2
#endif

     !> SGEHRD: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_sgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           real(sp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'SGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to zero
           do i = 1,ilo - 1
              tau(i) = zero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = zero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'SGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'SGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'SGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**t
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_slahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**t. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = one
                 call la_sgemm('NO TRANSPOSE','TRANSPOSE',ihi,ihi - i - ib + 1,ib,-one,work, &
                           ldwork,a(i + ib,i),lda,one,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_strmm('RIGHT','LOWER','TRANSPOSE','UNIT',i,ib - 1,one,a(i + 1,i) &
                           ,lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_saxpy(i,-one,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_slarfb('LEFT','TRANSPOSE','FORWARD','COLUMNWISE',ihi - i,n - i - ib + 1, &
                           ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work,ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_sgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_sgehrd
     !> DGEHRD: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_dgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           real(dp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'DGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to zero
           do i = 1,ilo - 1
              tau(i) = zero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = zero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'DGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'DGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'DGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**t
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_dlahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**t. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = one
                 call la_dgemm('NO TRANSPOSE','TRANSPOSE',ihi,ihi - i - ib + 1,ib,-one,work, &
                           ldwork,a(i + ib,i),lda,one,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_dtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',i,ib - 1,one,a(i + 1,i) &
                           ,lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_daxpy(i,-one,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_dlarfb('LEFT','TRANSPOSE','FORWARD','COLUMNWISE',ihi - i,n - i - ib + 1, &
                           ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work,ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_dgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_dgehrd
#ifdef LA_WITH_XDP
     !> XGEHRD: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_xgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           real(xdp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'XGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to zero
           do i = 1,ilo - 1
              tau(i) = zero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = zero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'XGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'XGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'XGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**t
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_xlahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**t. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = one
                 call la_xgemm('NO TRANSPOSE','TRANSPOSE',ihi,ihi - i - ib + 1,ib,-one,work, &
                           ldwork,a(i + ib,i),lda,one,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_xtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',i,ib - 1,one,a(i + 1,i) &
                           ,lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_xaxpy(i,-one,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_xlarfb('LEFT','TRANSPOSE','FORWARD','COLUMNWISE',ihi - i,n - i - ib + 1, &
                           ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work,ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_xgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_xgehrd
#endif
#ifdef LA_WITH_QP
     !> QGEHRD: reduces a real general matrix A to upper Hessenberg form H by
     !> an orthogonal similarity transformation:  Q**T * A * Q = H .

     pure subroutine la_qgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           real(qp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'QGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to zero
           do i = 1,ilo - 1
              tau(i) = zero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = zero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'QGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'QGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'QGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**t
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_qlahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**t. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = one
                 call la_qgemm('NO TRANSPOSE','TRANSPOSE',ihi,ihi - i - ib + 1,ib,-one,work, &
                           ldwork,a(i + ib,i),lda,one,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_qtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',i,ib - 1,one,a(i + 1,i) &
                           ,lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_qaxpy(i,-one,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_qlarfb('LEFT','TRANSPOSE','FORWARD','COLUMNWISE',ihi - i,n - i - ib + 1, &
                           ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work,ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_qgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_qgehrd
#endif

     !> CGEBAK: forms the right or left eigenvectors of a complex general
     !> matrix by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by CGEBAL.

     pure subroutine la_cgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(in) :: scale(*)
           complex(sp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(sp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('CGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_csscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_csscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_cswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_cswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_cgebak
     !> ZGEBAK: forms the right or left eigenvectors of a complex general
     !> matrix by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by ZGEBAL.

     pure subroutine la_zgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(in) :: scale(*)
           complex(dp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(dp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_zdscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_zdscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_zswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_zswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_zgebak
#ifdef LA_WITH_XDP
     !> YGEBAK: forms the right or left eigenvectors of a complex general
     !> matrix by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by YGEBAL.

     pure subroutine la_ygebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(in) :: scale(*)
           complex(xdp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(xdp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('YGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_yxscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_yxscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_yswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_yswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_ygebak
#endif
#ifdef LA_WITH_QP
     !> WGEBAK: forms the right or left eigenvectors of a complex general
     !> matrix by backward transformation on the computed eigenvectors of the
     !> balanced matrix output by WGEBAL.

     pure subroutine la_wgebak(job,side,n,ilo,ihi,scale,m,v,ldv,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job,side
           integer(ilp),intent(in) :: ihi,ilo,ldv,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(in) :: scale(*)
           complex(qp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: leftv,rightv
           integer(ilp) :: i,ii,k
           real(qp) :: s
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! decode and test the input parameters
           rightv = la_lsame(side,'R')
           leftv = la_lsame(side,'L')
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (.not. rightv .and. .not. leftv) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (m < 0) then
              info = -7
           else if (ldv < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WGEBAK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (m == 0) return
           if (la_lsame(job,'N')) return
           if (ilo == ihi) go to 30
           ! backward balance
           if (la_lsame(job,'S') .or. la_lsame(job,'B')) then
              if (rightv) then
                 do i = ilo,ihi
                    s = scale(i)
                    call la_wqscal(m,s,v(i,1),ldv)
                 end do
              end if
              if (leftv) then
                 do i = ilo,ihi
                    s = one/scale(i)
                    call la_wqscal(m,s,v(i,1),ldv)
                 end do
              end if
           end if
           ! backward permutation
           ! for  i = ilo-1 step -1 until 1,
                    ! ihi+1 step 1 until n do --
                    30 continue
           if (la_lsame(job,'P') .or. la_lsame(job,'B')) then
              if (rightv) then
                 loop_40: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_40
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_40
                    call la_wswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_40
              end if
              if (leftv) then
                 loop_50: do ii = 1,n
                    i = ii
                    if (i >= ilo .and. i <= ihi) cycle loop_50
                    if (i < ilo) i = ilo - ii
                    k = scale(i)
                    if (k == i) cycle loop_50
                    call la_wswap(m,v(i,1),ldv,v(k,1),ldv)
                 end do loop_50
              end if
           end if
           return
     end subroutine la_wgebak
#endif

     !> CGEBAL: balances a general complex matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_cgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: scale(*)
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sclfac = 2.0e+0_sp
           real(sp),parameter :: factor = 0.95e+0_sp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(sp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,min,real
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_cswap(l,a(1,j),1,a(1,m),1)
           call la_cswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (real(a(j,i),KIND=sp) /= zero .or. aimag(a(j,i)) /= zero) cycle &
                           loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (real(a(i,j),KIND=sp) /= zero .or. aimag(a(i,j)) /= zero) cycle &
                           loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_slamch('S')/la_slamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_scnrm2(l - k + 1,a(k,i),1)
              r = la_scnrm2(l - k + 1,a(i,k),lda)
              ica = la_icamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_icamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_sisnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('CGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_csscal(n - k + 1,g,a(i,k),lda)
              call la_csscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_cgebal
     !> ZGEBAL: balances a general complex matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_zgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: scale(*)
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sclfac = 2.0e+0_dp
           real(dp),parameter :: factor = 0.95e+0_dp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(dp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_zswap(l,a(1,j),1,a(1,m),1)
           call la_zswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (real(a(j,i),KIND=dp) /= zero .or. aimag(a(j,i)) /= zero) cycle &
                           loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (real(a(i,j),KIND=dp) /= zero .or. aimag(a(i,j)) /= zero) cycle &
                           loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_dlamch('S')/la_dlamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_dznrm2(l - k + 1,a(k,i),1)
              r = la_dznrm2(l - k + 1,a(i,k),lda)
              ica = la_izamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_izamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_disnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('ZGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_zdscal(n - k + 1,g,a(i,k),lda)
              call la_zdscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_zgebal
#ifdef LA_WITH_XDP
     !> YGEBAL: balances a general complex matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_ygebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: scale(*)
           complex(xdp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sclfac = 2.0e+0_xdp
           real(xdp),parameter :: factor = 0.95e+0_xdp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(xdp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('YGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_yswap(l,a(1,j),1,a(1,m),1)
           call la_yswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (real(a(j,i),KIND=xdp) /= zero .or. aimag(a(j,i)) /= zero) cycle &
                           loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (real(a(i,j),KIND=xdp) /= zero .or. aimag(a(i,j)) /= zero) cycle &
                           loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_xlamch('S')/la_xlamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_xynrm2(l - k + 1,a(k,i),1)
              r = la_xynrm2(l - k + 1,a(i,k),lda)
              ica = la_iyamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_iyamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_xisnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('YGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_yxscal(n - k + 1,g,a(i,k),lda)
              call la_yxscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_ygebal
#endif
#ifdef LA_WITH_QP
     !> WGEBAL: balances a general complex matrix A.  This involves, first,
     !> permuting A by a similarity transformation to isolate eigenvalues
     !> in the first 1 to ILO-1 and last IHI+1 to N elements on the
     !> diagonal; and second, applying a diagonal similarity transformation
     !> to rows and columns ILO to IHI to make the rows and columns as
     !> close in norm as possible.  Both steps are optional.
     !> Balancing may reduce the 1-norm of the matrix, and improve the
     !> accuracy of the computed eigenvalues and/or eigenvectors.

     pure subroutine la_wgebal(job,n,a,lda,ilo,ihi,scale,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: job
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: scale(*)
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sclfac = 2.0e+0_qp
           real(qp),parameter :: factor = 0.95e+0_qp

           ! Local Scalars
           logical(lk) :: noconv
           integer(ilp) :: i,ica,iexc,ira,j,k,l,m
           real(qp) :: c,ca,f,g,r,ra,s,sfmax1,sfmax2,sfmin1,sfmin2
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! test the input parameters
           info = 0
           if (.not. la_lsame(job,'N') .and. .not. la_lsame(job,'P') &
                     .and. .not. la_lsame(job,'S') .and. .not. la_lsame(job,'B')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WGEBAL',-info)
              return
           end if
           k = 1
           l = n
           if (n == 0) go to 210
           if (la_lsame(job,'N')) then
              do i = 1,n
                 scale(i) = one
              end do
              go to 210
           end if
           if (la_lsame(job,'S')) go to 120
           ! permutation to isolate eigenvalues if possible
           go to 50
           ! row and column exchange.
           20 continue
           scale(m) = j
           if (j == m) go to 30
           call la_wswap(l,a(1,j),1,a(1,m),1)
           call la_wswap(n - k + 1,a(j,k),lda,a(m,k),lda)
           30 continue
           go to(40,80) iexc
           ! search for rows isolating an eigenvalue and push them down.
           40 continue
           if (l == 1) go to 210
           l = l - 1
           50 continue
           loop_70: do j = l,1,-1
              loop_60: do i = 1,l
                 if (i == j) cycle loop_60
                 if (real(a(j,i),KIND=qp) /= zero .or. aimag(a(j,i)) /= zero) cycle &
                           loop_70
              end do loop_60
              m = l
              iexc = 1
              go to 20
           end do loop_70
           go to 90
           ! search for columns isolating an eigenvalue and push them left.
           80 continue
           k = k + 1
           90 continue
           loop_110: do j = k,l
              loop_100: do i = k,l
                 if (i == j) cycle loop_100
                 if (real(a(i,j),KIND=qp) /= zero .or. aimag(a(i,j)) /= zero) cycle &
                           loop_110
              end do loop_100
              m = k
              iexc = 2
              go to 20
           end do loop_110
           120 continue
           do i = k,l
              scale(i) = one
           end do
           if (la_lsame(job,'P')) go to 210
           ! balance the submatrix in rows k to l.
           ! iterative loop for norm reduction
           sfmin1 = la_qlamch('S')/la_qlamch('P')
           sfmax1 = one/sfmin1
           sfmin2 = sfmin1*sclfac
           sfmax2 = one/sfmin2
           140 continue
           noconv = .false.
           loop_200: do i = k,l
              c = la_qwnrm2(l - k + 1,a(k,i),1)
              r = la_qwnrm2(l - k + 1,a(i,k),lda)
              ica = la_iwamax(l,a(1,i),1)
              ca = abs(a(ica,i))
              ira = la_iwamax(n - k + 1,a(i,k),lda)
              ra = abs(a(i,ira + k - 1))
              ! guard against zero c or r due to underflow.
              if (c == zero .or. r == zero) cycle loop_200
              g = r/sclfac
              f = one
              s = c + r
              160 continue
              if (c >= g .or. max(f,c,ca) >= sfmax2 .or. min(r,g,ra) <= sfmin2) go to 170
                 if (la_qisnan(c + f + ca + r + g + ra)) then
                 ! exit if nan to avoid infinite loop
                 info = -3
                 call la_xerbla('WGEBAL',-info)
                 return
              end if
              f = f*sclfac
              c = c*sclfac
              ca = ca*sclfac
              r = r/sclfac
              g = g/sclfac
              ra = ra/sclfac
              go to 160
              170 continue
              g = c/sclfac
              180 continue
              if (g < r .or. max(r,ra) >= sfmax2 .or. min(f,c,g,ca) <= sfmin2) go to 190
              f = f/sclfac
              c = c/sclfac
              g = g/sclfac
              ca = ca/sclfac
              r = r*sclfac
              ra = ra*sclfac
              go to 180
              ! now balance.
              190 continue
              if ((c + r) >= factor*s) cycle loop_200
              if (f < one .and. scale(i) < one) then
                 if (f*scale(i) <= sfmin1) cycle loop_200
              end if
              if (f > one .and. scale(i) > one) then
                 if (scale(i) >= sfmax1/f) cycle loop_200
              end if
              g = one/f
              scale(i) = scale(i)*f
              noconv = .true.
              call la_wqscal(n - k + 1,g,a(i,k),lda)
              call la_wqscal(l,f,a(1,i),1)
           end do loop_200
           if (noconv) go to 140
           210 continue
           ilo = k
           ihi = l
           return
     end subroutine la_wgebal
#endif

     !> CGEHD2: reduces a complex general matrix A to upper Hessenberg form H
     !> by a unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_cgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              alpha = a(i + 1,i)
              call la_clarfg(ihi - i,alpha,a(min(i + 2,n),i),1,tau(i))
              a(i + 1,i) = cone
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_clarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i)**h to a(i+1:ihi,i+1:n) from the left
              call la_clarf('LEFT',ihi - i,n - i,a(i + 1,i),1,conjg(tau(i)),a(i + 1,i + &
                        1),lda,work)
              a(i + 1,i) = alpha
           end do
           return
     end subroutine la_cgehd2
     !> ZGEHD2: reduces a complex general matrix A to upper Hessenberg form H
     !> by a unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_zgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              alpha = a(i + 1,i)
              call la_zlarfg(ihi - i,alpha,a(min(i + 2,n),i),1,tau(i))
              a(i + 1,i) = cone
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_zlarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i)**h to a(i+1:ihi,i+1:n) from the left
              call la_zlarf('LEFT',ihi - i,n - i,a(i + 1,i),1,conjg(tau(i)),a(i + 1,i + &
                        1),lda,work)
              a(i + 1,i) = alpha
           end do
           return
     end subroutine la_zgehd2
#ifdef LA_WITH_XDP
     !> YGEHD2: reduces a complex general matrix A to upper Hessenberg form H
     !> by a unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_ygehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(xdp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('YGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              alpha = a(i + 1,i)
              call la_ylarfg(ihi - i,alpha,a(min(i + 2,n),i),1,tau(i))
              a(i + 1,i) = cone
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_ylarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i)**h to a(i+1:ihi,i+1:n) from the left
              call la_ylarf('LEFT',ihi - i,n - i,a(i + 1,i),1,conjg(tau(i)),a(i + 1,i + &
                        1),lda,work)
              a(i + 1,i) = alpha
           end do
           return
     end subroutine la_ygehd2
#endif
#ifdef LA_WITH_QP
     !> WGEHD2: reduces a complex general matrix A to upper Hessenberg form H
     !> by a unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_wgehd2(n,ilo,ihi,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WGEHD2',-info)
              return
           end if
           do i = ilo,ihi - 1
              ! compute elementary reflector h(i) to annihilate a(i+2:ihi,i)
              alpha = a(i + 1,i)
              call la_wlarfg(ihi - i,alpha,a(min(i + 2,n),i),1,tau(i))
              a(i + 1,i) = cone
              ! apply h(i) to a(1:ihi,i+1:ihi) from the right
              call la_wlarf('RIGHT',ihi,ihi - i,a(i + 1,i),1,tau(i),a(1,i + 1),lda, &
                        work)
              ! apply h(i)**h to a(i+1:ihi,i+1:n) from the left
              call la_wlarf('LEFT',ihi - i,n - i,a(i + 1,i),1,conjg(tau(i)),a(i + 1,i + &
                        1),lda,work)
              a(i + 1,i) = alpha
           end do
           return
     end subroutine la_wgehd2
#endif

     !> CLAHR2: reduces the first NB columns of A complex general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an unitary similarity transformation
     !> Q**H * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*v**H, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by CGEHRD.

     pure subroutine la_clahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(sp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**h
                 call la_clacgv(i - 1,a(k + i - 1,1),lda)
                 call la_cgemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,cone,a(k + 1,i),1)
                 call la_clacgv(i - 1,a(k + i - 1,1),lda)
                 ! apply i - v * t**h * v**h to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**h * b1
                 call la_ccopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_ctrmv('LOWER','CONJUGATE TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda, &
                           t(1,nb),1)
                 ! w := w + v2**h * b2
                 call la_cgemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a( &
                           k + i,i),1,cone,t(1,nb),1)
                 ! w := t**h * w
                 call la_ctrmv('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1, &
                           nb),1)
                 ! b2 := b2 - v2*w
                 call la_cgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-cone,a(k + i,1),lda,t(1,nb &
                           ),1,cone,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_ctrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_caxpy(i - 1,-cone,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_clarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = cone
              ! compute  y(k+1:n,i)
              call la_cgemv('NO TRANSPOSE',n - k,n - k - i + 1,cone,a(k + 1,i + 1),lda,a(k + i,i) &
                        ,1,czero,y(k + 1,i),1)
              call la_cgemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a(k + &
                        i,i),1,czero,t(1,i),1)
              call la_cgemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,t(1,i),1, &
                        cone,y(k + 1,i),1)
              call la_cscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_cscal(i - 1,-tau(i),t(1,i),1)
              call la_ctrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_clacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_ctrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,cone,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,cone,a(1, &
                      2 + nb),lda,a(k + 1 + nb,1),lda,cone,y,ldy)
           call la_ctrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,cone,t,ldt,y, &
                     ldy)
           return
     end subroutine la_clahr2
     !> ZLAHR2: reduces the first NB columns of A complex general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an unitary similarity transformation
     !> Q**H * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by ZGEHRD.

     pure subroutine la_zlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(dp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**h
                 call la_zlacgv(i - 1,a(k + i - 1,1),lda)
                 call la_zgemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,cone,a(k + 1,i),1)
                 call la_zlacgv(i - 1,a(k + i - 1,1),lda)
                 ! apply i - v * t**h * v**h to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**h * b1
                 call la_zcopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_ztrmv('LOWER','CONJUGATE TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda, &
                           t(1,nb),1)
                 ! w := w + v2**h * b2
                 call la_zgemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a( &
                           k + i,i),1,cone,t(1,nb),1)
                 ! w := t**h * w
                 call la_ztrmv('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1, &
                           nb),1)
                 ! b2 := b2 - v2*w
                 call la_zgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-cone,a(k + i,1),lda,t(1,nb &
                           ),1,cone,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_ztrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_zaxpy(i - 1,-cone,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_zlarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = cone
              ! compute  y(k+1:n,i)
              call la_zgemv('NO TRANSPOSE',n - k,n - k - i + 1,cone,a(k + 1,i + 1),lda,a(k + i,i) &
                        ,1,czero,y(k + 1,i),1)
              call la_zgemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a(k + &
                        i,i),1,czero,t(1,i),1)
              call la_zgemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,t(1,i),1, &
                        cone,y(k + 1,i),1)
              call la_zscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_zscal(i - 1,-tau(i),t(1,i),1)
              call la_ztrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_zlacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_ztrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,cone,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,cone,a(1, &
                      2 + nb),lda,a(k + 1 + nb,1),lda,cone,y,ldy)
           call la_ztrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,cone,t,ldt,y, &
                     ldy)
           return
     end subroutine la_zlahr2
#ifdef LA_WITH_XDP
     !> YLAHR2: reduces the first NB columns of A complex general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an unitary similarity transformation
     !> Q**H * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by YGEHRD.

     pure subroutine la_ylahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(xdp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**h
                 call la_ylacgv(i - 1,a(k + i - 1,1),lda)
                 call la_ygemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,cone,a(k + 1,i),1)
                 call la_ylacgv(i - 1,a(k + i - 1,1),lda)
                 ! apply i - v * t**h * v**h to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**h * b1
                 call la_ycopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_ytrmv('LOWER','CONJUGATE TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda, &
                           t(1,nb),1)
                 ! w := w + v2**h * b2
                 call la_ygemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a( &
                           k + i,i),1,cone,t(1,nb),1)
                 ! w := t**h * w
                 call la_ytrmv('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1, &
                           nb),1)
                 ! b2 := b2 - v2*w
                 call la_ygemv('NO TRANSPOSE',n - k - i + 1,i - 1,-cone,a(k + i,1),lda,t(1,nb &
                           ),1,cone,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_ytrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_yaxpy(i - 1,-cone,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_ylarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = cone
              ! compute  y(k+1:n,i)
              call la_ygemv('NO TRANSPOSE',n - k,n - k - i + 1,cone,a(k + 1,i + 1),lda,a(k + i,i) &
                        ,1,czero,y(k + 1,i),1)
              call la_ygemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a(k + &
                        i,i),1,czero,t(1,i),1)
              call la_ygemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,t(1,i),1, &
                        cone,y(k + 1,i),1)
              call la_yscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_yscal(i - 1,-tau(i),t(1,i),1)
              call la_ytrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_ylacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_ytrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,cone,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,cone,a(1, &
                      2 + nb),lda,a(k + 1 + nb,1),lda,cone,y,ldy)
           call la_ytrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,cone,t,ldt,y, &
                     ldy)
           return
     end subroutine la_ylahr2
#endif
#ifdef LA_WITH_QP
     !> WLAHR2: reduces the first NB columns of A complex general n-BY-(n-k+1)
     !> matrix A so that elements below the k-th subdiagonal are zero. The
     !> reduction is performed by an unitary similarity transformation
     !> Q**H * A * Q. The routine returns the matrices V and T which determine
     !> Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T.
     !> This is an auxiliary routine called by WGEHRD.

     pure subroutine la_wlahr2(n,k,nb,a,lda,tau,t,ldt,y,ldy)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: k,lda,ldt,ldy,n,nb
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: t(ldt,nb),tau(nb),y(ldy,nb)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(qp) :: ei
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           ! quick return if possible
           if (n <= 1) return
           loop_10: do i = 1,nb
              if (i > 1) then
                 ! update a(k+1:n,i)
                 ! update i-th column of a - y * v**h
                 call la_wlacgv(i - 1,a(k + i - 1,1),lda)
                 call la_wgemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,a(k + i - 1,1), &
                           lda,cone,a(k + 1,i),1)
                 call la_wlacgv(i - 1,a(k + i - 1,1),lda)
                 ! apply i - v * t**h * v**h to this column (call it b) from the
                 ! left, using the last column of t as workspace
                 ! let  v = ( v1 )   and   b = ( b1 )   (first i-1 rows)
                          ! ( v2 )             ( b2 )
                 ! where v1 is unit lower triangular
                 ! w := v1**h * b1
                 call la_wcopy(i - 1,a(k + 1,i),1,t(1,nb),1)
                 call la_wtrmv('LOWER','CONJUGATE TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda, &
                           t(1,nb),1)
                 ! w := w + v2**h * b2
                 call la_wgemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a( &
                           k + i,i),1,cone,t(1,nb),1)
                 ! w := t**h * w
                 call la_wtrmv('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1, &
                           nb),1)
                 ! b2 := b2 - v2*w
                 call la_wgemv('NO TRANSPOSE',n - k - i + 1,i - 1,-cone,a(k + i,1),lda,t(1,nb &
                           ),1,cone,a(k + i,i),1)
                 ! b1 := b1 - v1*w
                 call la_wtrmv('LOWER','NO TRANSPOSE','UNIT',i - 1,a(k + 1,1),lda,t(1, &
                           nb),1)
                 call la_waxpy(i - 1,-cone,t(1,nb),1,a(k + 1,i),1)
                 a(k + i - 1,i - 1) = ei
              end if
              ! generate the elementary reflector h(i) to annihilate
              ! a(k+i+1:n,i)
              call la_wlarfg(n - k - i + 1,a(k + i,i),a(min(k + i + 1,n),i),1,tau(i))

              ei = a(k + i,i)
              a(k + i,i) = cone
              ! compute  y(k+1:n,i)
              call la_wgemv('NO TRANSPOSE',n - k,n - k - i + 1,cone,a(k + 1,i + 1),lda,a(k + i,i) &
                        ,1,czero,y(k + 1,i),1)
              call la_wgemv('CONJUGATE TRANSPOSE',n - k - i + 1,i - 1,cone,a(k + i,1),lda,a(k + &
                        i,i),1,czero,t(1,i),1)
              call la_wgemv('NO TRANSPOSE',n - k,i - 1,-cone,y(k + 1,1),ldy,t(1,i),1, &
                        cone,y(k + 1,i),1)
              call la_wscal(n - k,tau(i),y(k + 1,i),1)
              ! compute t(1:i,i)
              call la_wscal(i - 1,-tau(i),t(1,i),1)
              call la_wtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i),1)

              t(i,i) = tau(i)
           end do loop_10
           a(k + nb,nb) = ei
           ! compute y(1:k,1:nb)
           call la_wlacpy('ALL',k,nb,a(1,2),lda,y,ldy)
           call la_wtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',k,nb,cone,a(k + 1,1), &
                     lda,y,ldy)
           if (n > k + nb) call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',k,nb,n - k - nb,cone,a(1, &
                      2 + nb),lda,a(k + 1 + nb,1),lda,cone,y,ldy)
           call la_wtrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',k,nb,cone,t,ldt,y, &
                     ldy)
           return
     end subroutine la_wlahr2
#endif

     !> CUNGHR: generates a complex unitary matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> CGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_cunghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'CUNGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CUNGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors cone
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = czero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = czero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_cungqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_cunghr
     !> ZUNGHR: generates a complex unitary matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> ZGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_zunghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'ZUNGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZUNGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors cone
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = czero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = czero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_zungqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_zunghr
#ifdef LA_WITH_XDP
     !> YUNGHR: generates a complex unitary matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> YGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_yunghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'YUNGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YUNGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors cone
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = czero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = czero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_yungqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_yunghr
#endif
#ifdef LA_WITH_QP
     !> WUNGHR: generates a complex unitary matrix Q which is defined as the
     !> product of IHI-ILO elementary reflectors of order N, as returned by
     !> WGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_wunghr(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,lwkopt,nb,nh
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,nh) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              nb = la_ilaenv(1,'WUNGQR',' ',nh,nh,nh,-1)
              lwkopt = max(1,nh)*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WUNGHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              work(1) = 1
              return
           end if
           ! shift the vectors which define the elementary reflectors cone
           ! column to the right, and set the first ilo and the last n-ihi
           ! rows and columns to those of the unit matrix
           do j = ihi,ilo + 1,-1
              do i = 1,j - 1
                 a(i,j) = czero
              end do
              do i = j + 1,ihi
                 a(i,j) = a(i,j - 1)
              end do
              do i = ihi + 1,n
                 a(i,j) = czero
              end do
           end do
           do j = 1,ilo
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           do j = ihi + 1,n
              do i = 1,n
                 a(i,j) = czero
              end do
              a(j,j) = cone
           end do
           if (nh > 0) then
              ! generate q(ilo+1:ihi,ilo+1:ihi)
              call la_wungqr(nh,nh,nh,a(ilo + 1,ilo + 1),lda,tau(ilo),work,lwork, &
                        iinfo)
           end if
           work(1) = lwkopt
           return
     end subroutine la_wunghr
#endif

     !> CUNMHR: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by CGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_cunmhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'CUNMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'CUNMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CUNMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_cunmqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_cunmhr
     !> ZUNMHR: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by ZGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_zunmhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'ZUNMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'ZUNMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZUNMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_zunmqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_zunmhr
#ifdef LA_WITH_XDP
     !> YUNMHR: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by YGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_yunmhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'YUNMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'YUNMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YUNMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_yunmqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_yunmhr
#endif
#ifdef LA_WITH_QP
     !> WUNMHR: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix of order nq, with nq = m if
     !> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
     !> IHI-ILO elementary reflectors, as returned by WGEHRD:
     !> Q = H(ilo) H(ilo+1) . . . H(ihi-1).

     pure subroutine la_wunmhr(side,trans,m,n,ilo,ihi,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: ihi,ilo,lda,ldc,lwork,m,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,lquery
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,nh,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nh = ihi - ilo
           left = la_lsame(side,'L')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (ilo < 1 .or. ilo > max(1,nq)) then
              info = -5
           else if (ihi < min(ilo,nq) .or. ihi > nq) then
              info = -6
           else if (lda < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (left) then
                 nb = la_ilaenv(1,'WUNMQR',side//trans,nh,n,nh,-1)
              else
                 nb = la_ilaenv(1,'WUNMQR',side//trans,m,nh,nh,-1)
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WUNMHR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. nh == 0) then
              work(1) = 1
              return
           end if
           if (left) then
              mi = nh
              ni = n
              i1 = ilo + 1
              i2 = 1
           else
              mi = m
              ni = nh
              i1 = 1
              i2 = ilo + 1
           end if
           call la_wunmqr(side,trans,mi,ni,nh,a(ilo + 1,ilo),lda,tau(ilo),c(i1, &
                     i2),ldc,work,lwork,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_wunmhr
#endif

     !> CGEHRD: reduces a complex general matrix A to upper Hessenberg form H by
     !> an unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_cgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           complex(sp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'CGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to czero
           do i = 1,ilo - 1
              tau(i) = czero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = czero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'CGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'CGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'CGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**h
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_clahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**h. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = cone
                 call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',ihi,ihi - i - ib + 1,ib,- &
                           cone,work,ldwork,a(i + ib,i),lda,cone,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_ctrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',i,ib - 1,cone, &
                           a(i + 1,i),lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_caxpy(i,-cone,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_clarfb('LEFT','CONJUGATE TRANSPOSE','FORWARD','COLUMNWISE',ihi - i, &
                 n - i - ib + 1,ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work, &
                           ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_cgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_cgehrd
     !> ZGEHRD: reduces a complex general matrix A to upper Hessenberg form H by
     !> an unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_zgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           complex(dp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'ZGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to czero
           do i = 1,ilo - 1
              tau(i) = czero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = czero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'ZGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'ZGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'ZGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**h
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_zlahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**h. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = cone
                 call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',ihi,ihi - i - ib + 1,ib,- &
                           cone,work,ldwork,a(i + ib,i),lda,cone,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_ztrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',i,ib - 1,cone, &
                           a(i + 1,i),lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_zaxpy(i,-cone,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_zlarfb('LEFT','CONJUGATE TRANSPOSE','FORWARD','COLUMNWISE',ihi - i, &
                 n - i - ib + 1,ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work, &
                           ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_zgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_zgehrd
#ifdef LA_WITH_XDP
     !> YGEHRD: reduces a complex general matrix A to upper Hessenberg form H by
     !> an unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_ygehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           complex(xdp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'YGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to czero
           do i = 1,ilo - 1
              tau(i) = czero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = czero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'YGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'YGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'YGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**h
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_ylahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**h. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = cone
                 call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',ihi,ihi - i - ib + 1,ib,- &
                           cone,work,ldwork,a(i + ib,i),lda,cone,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_ytrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',i,ib - 1,cone, &
                           a(i + 1,i),lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_yaxpy(i,-cone,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_ylarfb('LEFT','CONJUGATE TRANSPOSE','FORWARD','COLUMNWISE',ihi - i, &
                 n - i - ib + 1,ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work, &
                           ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_ygehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_ygehrd
#endif
#ifdef LA_WITH_QP
     !> WGEHRD: reduces a complex general matrix A to upper Hessenberg form H by
     !> an unitary similarity transformation:  Q**H * A * Q = H .

     pure subroutine la_wgehrd(n,ilo,ihi,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,lda,lwork,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iwt,j,ldwork,lwkopt,nb,nbmin,nh,nx
           complex(qp) :: ei
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -2
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'WGEHRD',' ',n,ilo,ihi,-1))
              lwkopt = n*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WGEHRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! set elements 1:ilo-1 and ihi:n-1 of tau to czero
           do i = 1,ilo - 1
              tau(i) = czero
           end do
           do i = max(1,ihi),n - 1
              tau(i) = czero
           end do
           ! quick return if possible
           nh = ihi - ilo + 1
           if (nh <= 1) then
              work(1) = 1
              return
           end if
           ! determine the block size
           nb = min(nbmax,la_ilaenv(1,'WGEHRD',' ',n,ilo,ihi,-1))
           nbmin = 2
           if (nb > 1 .and. nb < nh) then
              ! determine when to cross over from blocked to unblocked code
              ! (last block is always handled by unblocked code)
              nx = max(nb,la_ilaenv(3,'WGEHRD',' ',n,ilo,ihi,-1))
              if (nx < nh) then
                 ! determine if workspace is large enough for blocked code
                 if (lwork < n*nb + tsize) then
                    ! not enough workspace to use optimal nb:  determine the
                    ! minimum value of nb, and reduce nb or force use of
                    ! unblocked code
                    nbmin = max(2,la_ilaenv(2,'WGEHRD',' ',n,ilo,ihi,-1))
                    if (lwork >= (n*nbmin + tsize)) then
                       nb = (lwork - tsize)/n
                    else
                       nb = 1
                    end if
                 end if
              end if
           end if
           ldwork = n
           if (nb < nbmin .or. nb >= nh) then
              ! use unblocked code below
              i = ilo
           else
              ! use blocked code
              iwt = 1 + n*nb
              do i = ilo,ihi - 1 - nx,nb
                 ib = min(nb,ihi - i)
                 ! reduce columns i:i+ib-1 to hessenberg form, returning the
                 ! matrices v and t of the block reflector h = i - v*t*v**h
                 ! which performs the reduction, and also the matrix y = a*v*t
                 call la_wlahr2(ihi,i,ib,a(1,i),lda,tau(i),work(iwt),ldt,work, &
                           ldwork)
                 ! apply the block reflector h to a(1:ihi,i+ib:ihi) from the
                 ! right, computing  a := a - y * v**h. v(i+ib,ib-1) must be set
                 ! to 1
                 ei = a(i + ib,i + ib - 1)
                 a(i + ib,i + ib - 1) = cone
                 call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',ihi,ihi - i - ib + 1,ib,- &
                           cone,work,ldwork,a(i + ib,i),lda,cone,a(1,i + ib),lda)
                 a(i + ib,i + ib - 1) = ei
                 ! apply the block reflector h to a(1:i,i+1:i+ib-1) from the
                 ! right
                 call la_wtrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',i,ib - 1,cone, &
                           a(i + 1,i),lda,work,ldwork)
                 do j = 0,ib - 2
                    call la_waxpy(i,-cone,work(ldwork*j + 1),1,a(1,i + j + 1),1)
                 end do
                 ! apply the block reflector h to a(i+1:ihi,i+ib:n) from the
                 ! left
                 call la_wlarfb('LEFT','CONJUGATE TRANSPOSE','FORWARD','COLUMNWISE',ihi - i, &
                 n - i - ib + 1,ib,a(i + 1,i),lda,work(iwt),ldt,a(i + 1,i + ib),lda,work, &
                           ldwork)
              end do
           end if
           ! use unblocked code to reduce the rest of the matrix
           call la_wgehd2(n,i,ihi,a,lda,tau,work,iinfo)
           work(1) = lwkopt
           return
     end subroutine la_wgehrd
#endif

end module la_lapack_eigv_gen_hess
