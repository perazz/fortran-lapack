!> Least-squares drivers: QR, complete orthogonal, SVD and divide-and-conquer solutions
module la_lapack_lsq
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_l2
     use la_lapack_blas_like_mnorm
     use la_lapack_lsq_aux
     use la_lapack_orthogonal_factors_ql
     use la_lapack_orthogonal_factors_qr
     use la_lapack_orthogonal_factors_rz
     use la_lapack_solve_tri_comp
     use la_lapack_svd_bidiag_qr
     use la_lapack_svd_comp
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sgels
     public :: la_sgelsy
     public :: la_sgetsls
     public :: la_sgelsd
     public :: la_sgelss
     public :: la_dgels
     public :: la_dgelsy
     public :: la_dgetsls
     public :: la_dgelsd
     public :: la_dgelss
#ifdef LA_WITH_XDP
     public :: la_xgels
     public :: la_xgelsy
     public :: la_xgetsls
     public :: la_xgelsd
     public :: la_xgelss
#endif
#ifdef LA_WITH_QP
     public :: la_qgels
     public :: la_qgelsy
     public :: la_qgetsls
     public :: la_qgelsd
     public :: la_qgelss
#endif
     public :: la_cgels
     public :: la_cgelsd
     public :: la_cgelss
     public :: la_cgelsy
     public :: la_cgetsls
     public :: la_zgels
     public :: la_zgelsd
     public :: la_zgelss
     public :: la_zgelsy
     public :: la_zgetsls
#ifdef LA_WITH_XDP
     public :: la_ygels
     public :: la_ygelsd
     public :: la_ygelss
     public :: la_ygelsy
     public :: la_ygetsls
#endif
#ifdef LA_WITH_QP
     public :: la_wgels
     public :: la_wgelsd
     public :: la_wgelss
     public :: la_wgelsy
     public :: la_wgetsls
#endif

     contains

     !> SGELS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, or its transpose, using a QR or LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_sgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(sp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(sp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'SGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'SORMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'SORMQR','LT',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'SGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'SORMLQ','LT',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'SORMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=sp)
           end if
           if (info /= 0) then
              call la_xerbla('SGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_slaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_slamch('S')/la_slamch('P')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_slange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_slaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_slange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_sgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
                 call la_sormqr('LEFT','TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_strtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_strtrs('UPPER','TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_sormqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_sgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_strtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_sormlq('LEFT','TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_sormlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_strtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_slascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_slascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_slascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_slascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=sp)
           return
     end subroutine la_sgels
     !> DGELS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, or its transpose, using a QR or LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(dp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(dp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'DGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'DORMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'DORMQR','LT',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'DGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'DORMLQ','LT',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'DORMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=dp)
           end if
           if (info /= 0) then
              call la_xerbla('DGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_dlaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_dlamch('S')/la_dlamch('P')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_dlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_dlange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_dgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
                 call la_dormqr('LEFT','TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_dtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_dtrtrs('UPPER','TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_dormqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_dgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_dtrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_dormlq('LEFT','TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_dormlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_dtrtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_dlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_dlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_dlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_dlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=dp)
           return
     end subroutine la_dgels
#ifdef LA_WITH_XDP
     !> XGELS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, or its transpose, using a QR or LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_xgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(xdp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(xdp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'XGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'XORMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'XORMQR','LT',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'XGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'XORMLQ','LT',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'XORMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=xdp)
           end if
           if (info /= 0) then
              call la_xerbla('XGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_xlaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_xlamch('S')/la_xlamch('P')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_xlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_xlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_xlange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_xgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
                 call la_xormqr('LEFT','TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_xtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_xtrtrs('UPPER','TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_xormqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_xgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_xtrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_xormlq('LEFT','TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_xormlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_xtrtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_xlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_xlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_xlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_xlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=xdp)
           return
     end subroutine la_xgels
#endif
#ifdef LA_WITH_QP
     !> QGELS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, or its transpose, using a QR or LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_qgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(qp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(qp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'QGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'QORMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'QORMQR','LT',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'QGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'QORMLQ','LT',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'QORMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=qp)
           end if
           if (info /= 0) then
              call la_xerbla('QGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_qlaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_qlamch('S')/la_qlamch('P')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_qlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_qlange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_qgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
                 call la_qormqr('LEFT','TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_qtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_qtrtrs('UPPER','TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_qormqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_qgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_qtrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_qormlq('LEFT','TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                           work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_qormlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_qtrtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_qlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_qlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_qlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_qlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=qp)
           return
     end subroutine la_qgels
#endif

     !> SGELSY: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by orthogonal transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**T [ inv(T11)*Q1**T*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.

     subroutine la_sgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,info)
        use la_constants_sp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(sp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkmin,lwkopt,mn,nb,nb1,nb2, &
                     nb3,nb4
           real(sp) :: anrm,bignum,bnrm,c1,c2,s1,s2,smax,smaxpr,smin,sminpr,smlnum, &
                     wsize
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           end if
           ! figure out optimal block size
           if (info == 0) then
              if (mn == 0 .or. nrhs == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'SGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'SGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'SORMQR',' ',m,n,nrhs,-1)
                 nb4 = la_ilaenv(1,'SORMRQ',' ',m,n,nrhs,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = mn + max(2*mn,n + 1,mn + nrhs)
                 lwkopt = max(lwkmin,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (mn == 0 .or. nrhs == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_slamch('S')/la_slamch('P')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_slange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_slaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_slange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_sgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,info)

           wsize = mn + work(mn + 1)
           ! workspace: mn+2*n+nb*(n+1).
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = one
           work(ismax) = one
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_slaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_slaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_slaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_stzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
           call la_sormqr('LEFT','TRANSPOSE',m,nrhs,mn,a,lda,work(1),b,ldb,work( &
                     2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + work(2*mn + 1))
           ! workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_strsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,one,a,lda, &
                      b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = zero
              end do
           end do
           ! b(1:n,1:nrhs) := y**t * b(1:n,1:nrhs)
           if (rank < n) then
              call la_sormrz('LEFT','TRANSPOSE',n,nrhs,rank,n - rank,a,lda,work(mn + 1), &
                         b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_scopy(n,work(1),1,b(1,j),1)
           end do
           ! workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_slascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_slascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_slascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_slascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_slascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_slascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = lwkopt
           return
     end subroutine la_sgelsy
     !> DGELSY: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by orthogonal transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**T [ inv(T11)*Q1**T*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.

     subroutine la_dgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,info)
        use la_constants_dp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(dp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkmin,lwkopt,mn,nb,nb1,nb2, &
                     nb3,nb4
           real(dp) :: anrm,bignum,bnrm,c1,c2,s1,s2,smax,smaxpr,smin,sminpr,smlnum, &
                     wsize
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           end if
           ! figure out optimal block size
           if (info == 0) then
              if (mn == 0 .or. nrhs == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'DGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'DGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'DORMQR',' ',m,n,nrhs,-1)
                 nb4 = la_ilaenv(1,'DORMRQ',' ',m,n,nrhs,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = mn + max(2*mn,n + 1,mn + nrhs)
                 lwkopt = max(lwkmin,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (mn == 0 .or. nrhs == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_dlamch('S')/la_dlamch('P')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_dlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_dlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_dlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_dgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,info)

           wsize = mn + work(mn + 1)
           ! workspace: mn+2*n+nb*(n+1).
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = one
           work(ismax) = one
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_dlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_dlaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_dlaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_dtzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
           call la_dormqr('LEFT','TRANSPOSE',m,nrhs,mn,a,lda,work(1),b,ldb,work( &
                     2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + work(2*mn + 1))
           ! workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_dtrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,one,a,lda, &
                      b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = zero
              end do
           end do
           ! b(1:n,1:nrhs) := y**t * b(1:n,1:nrhs)
           if (rank < n) then
              call la_dormrz('LEFT','TRANSPOSE',n,nrhs,rank,n - rank,a,lda,work(mn + 1), &
                         b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_dcopy(n,work(1),1,b(1,j),1)
           end do
           ! workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_dlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_dlascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_dlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_dlascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_dlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_dlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = lwkopt
           return
     end subroutine la_dgelsy
#ifdef LA_WITH_XDP
     !> XGELSY: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by orthogonal transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**T [ inv(T11)*Q1**T*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.

     subroutine la_xgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,info)
        use la_constants_xdp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(xdp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkmin,lwkopt,mn,nb,nb1,nb2, &
                     nb3,nb4
           real(xdp) :: anrm,bignum,bnrm,c1,c2,s1,s2,smax,smaxpr,smin,sminpr,smlnum, &
                     wsize
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           end if
           ! figure out optimal block size
           if (info == 0) then
              if (mn == 0 .or. nrhs == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'XGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'XGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'XORMQR',' ',m,n,nrhs,-1)
                 nb4 = la_ilaenv(1,'XORMRQ',' ',m,n,nrhs,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = mn + max(2*mn,n + 1,mn + nrhs)
                 lwkopt = max(lwkmin,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (mn == 0 .or. nrhs == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_xlamch('S')/la_xlamch('P')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_xlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_xlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_xlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_xgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,info)

           wsize = mn + work(mn + 1)
           ! workspace: mn+2*n+nb*(n+1).
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = one
           work(ismax) = one
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_xlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_xlaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_xlaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_xtzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
           call la_xormqr('LEFT','TRANSPOSE',m,nrhs,mn,a,lda,work(1),b,ldb,work( &
                     2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + work(2*mn + 1))
           ! workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_xtrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,one,a,lda, &
                      b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = zero
              end do
           end do
           ! b(1:n,1:nrhs) := y**t * b(1:n,1:nrhs)
           if (rank < n) then
              call la_xormrz('LEFT','TRANSPOSE',n,nrhs,rank,n - rank,a,lda,work(mn + 1), &
                         b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_xcopy(n,work(1),1,b(1,j),1)
           end do
           ! workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_xlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_xlascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_xlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_xlascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_xlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_xlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = lwkopt
           return
     end subroutine la_xgelsy
#endif
#ifdef LA_WITH_QP
     !> QGELSY: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by orthogonal transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**T [ inv(T11)*Q1**T*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.

     subroutine la_qgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,info)
        use la_constants_qp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(qp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkmin,lwkopt,mn,nb,nb1,nb2, &
                     nb3,nb4
           real(qp) :: anrm,bignum,bnrm,c1,c2,s1,s2,smax,smaxpr,smin,sminpr,smlnum, &
                     wsize
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           end if
           ! figure out optimal block size
           if (info == 0) then
              if (mn == 0 .or. nrhs == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'QGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'QGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'QORMQR',' ',m,n,nrhs,-1)
                 nb4 = la_ilaenv(1,'QORMRQ',' ',m,n,nrhs,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = mn + max(2*mn,n + 1,mn + nrhs)
                 lwkopt = max(lwkmin,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (mn == 0 .or. nrhs == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_qlamch('S')/la_qlamch('P')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_qlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_qlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_qlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_qgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,info)

           wsize = mn + work(mn + 1)
           ! workspace: mn+2*n+nb*(n+1).
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = one
           work(ismax) = one
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_qlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_qlaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_qlaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_qtzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
           call la_qormqr('LEFT','TRANSPOSE',m,nrhs,mn,a,lda,work(1),b,ldb,work( &
                     2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + work(2*mn + 1))
           ! workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_qtrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,one,a,lda, &
                      b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = zero
              end do
           end do
           ! b(1:n,1:nrhs) := y**t * b(1:n,1:nrhs)
           if (rank < n) then
              call la_qormrz('LEFT','TRANSPOSE',n,nrhs,rank,n - rank,a,lda,work(mn + 1), &
                         b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_qcopy(n,work(1),1,b(1,j),1)
           end do
           ! workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_qlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_qlascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_qlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_qlascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_qlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_qlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = lwkopt
           return
     end subroutine la_qgelsy
#endif

     !> SGETSLS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_sgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(sp) :: anrm,bignum,bnrm,smlnum,tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'T')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_sgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_sgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_sgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_sgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_sgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_sgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_sgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_sgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=sp)
           end if
           if (info /= 0) then
             call la_xerbla('SGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=sp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_slaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_slamch('S')/la_slamch('P')
            bignum = one/smlnum
            call la_slabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_slange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_slaset('F',maxmn,nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_slange('M',brow,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_sgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_sgemqr('L','T',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_strtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_strtrs('U','T','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_sgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_sgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_strtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_sgemlq('L','T',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_sgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_strtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_slascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_slascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_slascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_slascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=sp)
           return
     end subroutine la_sgetsls
     !> DGETSLS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_dgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(dp) :: anrm,bignum,bnrm,smlnum,tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'T')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_dgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_dgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_dgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_dgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_dgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_dgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_dgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_dgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=dp)
           end if
           if (info /= 0) then
             call la_xerbla('DGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=dp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_dlaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_dlamch('S')/la_dlamch('P')
            bignum = one/smlnum
            call la_dlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_dlaset('F',maxmn,nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_dlange('M',brow,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_dgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_dgemqr('L','T',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_dtrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_dtrtrs('U','T','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_dgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_dgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_dtrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_dgemlq('L','T',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_dgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_dtrtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_dlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_dlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_dlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_dlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=dp)
           return
     end subroutine la_dgetsls
#ifdef LA_WITH_XDP
     !> XGETSLS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_xgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(xdp) :: anrm,bignum,bnrm,smlnum,tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'T')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_xgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_xgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_xgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_xgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_xgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_xgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_xgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_xgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=xdp)
           end if
           if (info /= 0) then
             call la_xerbla('XGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=xdp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_xlaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_xlamch('S')/la_xlamch('P')
            bignum = one/smlnum
            call la_xlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_xlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_xlaset('F',maxmn,nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_xlange('M',brow,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_xgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_xgemqr('L','T',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_xtrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_xtrtrs('U','T','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_xgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_xgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_xtrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_xgemlq('L','T',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_xgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_xtrtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_xlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_xlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_xlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_xlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=xdp)
           return
     end subroutine la_xgetsls
#endif
#ifdef LA_WITH_QP
     !> QGETSLS: solves overdetermined or underdetermined real linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'T' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_qgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(qp) :: anrm,bignum,bnrm,smlnum,tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'T')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'T'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_qgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_qgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_qgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_qgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_qgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_qgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_qgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_qgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=qp)
           end if
           if (info /= 0) then
             call la_xerbla('QGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=qp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_qlaset('FULL',max(m,n),nrhs,zero,zero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_qlamch('S')/la_qlamch('P')
            bignum = one/smlnum
            call la_qlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_qlaset('F',maxmn,nrhs,zero,zero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_qlange('M',brow,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_qgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_qgemqr('L','T',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_qtrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_qtrtrs('U','T','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_qgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_qgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_qtrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = zero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_qgemlq('L','T',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_qgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_qtrtrs('LOWER','TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_qlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_qlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_qlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_qlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=qp)
           return
     end subroutine la_qgetsls
#endif

     !> SGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_sgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,iwork, &
               info)
        use la_constants_sp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(sp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,maxmn, &
                     maxwrk,minmn,minwrk,mm,mnthr,nlvl,nwork,smlsiz,wlalsd
           real(sp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: int,log,max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              liwork = 1
              if (minmn > 0) then
                 smlsiz = la_ilaenv(9,'SGELSD',' ',0,0,0,0)
                 mnthr = la_ilaenv(6,'SGELSD',' ',m,n,nrhs,-1)
                 nlvl = max(int(log(real(minmn,KIND=sp)/real(smlsiz + 1,KIND=sp))/log( &
                           two),KIND=ilp) + 1,0)
              liwork = 3*minmn*nlvl + 11*minmn
              mm = m
              if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns.
                 mm = n
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'SGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'SORMQR','LT',m,nrhs,n,- &
                              1))
              end if
              if (m >= n) then
                 ! path 1 - overdetermined or exactly determined.
                    maxwrk = max(maxwrk,3*n + (mm + n)*la_ilaenv(1,'SGEBRD',' ',mm,n, &
                              -1,-1))
                    maxwrk = max(maxwrk,3*n + nrhs*la_ilaenv(1,'SORMBR','QLT',mm,nrhs, &
                              n,-1))
                    maxwrk = max(maxwrk,3*n + (n - 1)*la_ilaenv(1,'SORMBR','PLN',n, &
                              nrhs,n,-1))
                 wlalsd = 9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + (smlsiz + 1)**2
                 maxwrk = max(maxwrk,3*n + wlalsd)
                 minwrk = max(3*n + mm,3*n + nrhs,3*n + wlalsd)
              end if
              if (n > m) then
                 wlalsd = 9*m + 2*m*smlsiz + 8*m*nlvl + m*nrhs + (smlsiz + 1)**2
                 if (n >= mnthr) then
                    ! path 2a - underdetermined, with many more columns
                    ! than rows.
                    maxwrk = m + m*la_ilaenv(1,'SGELQF',' ',m,n,-1,-1)
                       maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'SGEBRD',' ',m,m, &
                                  -1,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'SORMBR','QLT',m, &
                              nrhs,m,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'SORMBR', &
                                 'PLN',m,nrhs,m,-1))
                    if (nrhs > 1) then
                       maxwrk = max(maxwrk,m*m + m + m*nrhs)
                    else
                       maxwrk = max(maxwrk,m*m + 2*m)
                    end if
                       maxwrk = max(maxwrk,m + nrhs*la_ilaenv(1,'SORMLQ','LT',n,nrhs,m, &
                                  -1))
                    maxwrk = max(maxwrk,m*m + 4*m + wlalsd)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                    maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                 else
                    ! path 2 - remaining underdetermined cases.
                    maxwrk = 3*m + (n + m)*la_ilaenv(1,'SGEBRD',' ',m,n,-1,-1)

                       maxwrk = max(maxwrk,3*m + nrhs*la_ilaenv(1,'SORMBR','QLT',m,nrhs, &
                                  n,-1))
                       maxwrk = max(maxwrk,3*m + m*la_ilaenv(1,'SORMBR','PLN',n,nrhs,m, &
                              -1))
                    maxwrk = max(maxwrk,3*m + wlalsd)
                 end if
                 minwrk = max(3*m + nrhs,3*m + m,3*m + wlalsd)
                 end if
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGELSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_slamch('P')
           sfmin = la_slamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_slange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_slascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_slascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_slaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_slaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_slange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_slascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_slascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure certain entries of b are zero.
           if (m < n) call la_slaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_sgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_sormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_slaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
                 end if
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              nwork = itaup + n
              ! bidiagonalize r in a.
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_sgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_sormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_slalsd('U',smlsiz,n,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_sormbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m,wlalsd)) &
                     then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs,4*m + m*lda + &
                        wlalsd)) ldwork = lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_sgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_slacpy('L',m,m,a,lda,work(il),ldwork)
              call la_slaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize l in work(il).
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_sgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_sormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_slalsd('U',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_sormbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_slaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_sormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize a.
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_sgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_sormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_slalsd('L',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_sormbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_slascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_slascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_slascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_slascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           return
     end subroutine la_sgelsd
     !> DGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_dgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,iwork, &
               info)
        use la_constants_dp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(dp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,maxmn, &
                     maxwrk,minmn,minwrk,mm,mnthr,nlvl,nwork,smlsiz,wlalsd
           real(dp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: real,int,log,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           mnthr = la_ilaenv(6,'DGELSD',' ',m,n,nrhs,-1)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           smlsiz = la_ilaenv(9,'DGELSD',' ',0,0,0,0)
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           minwrk = 1
           liwork = 1
           minmn = max(1,minmn)
           nlvl = max(int(log(real(minmn,KIND=dp)/real(smlsiz + 1,KIND=dp))/log(two), &
                     KIND=ilp) + 1,0)
           if (info == 0) then
              maxwrk = 0
              liwork = 3*minmn*nlvl + 11*minmn
              mm = m
              if (m >= n .and. m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'DGEQRF',' ',m,n,-1,-1))

                 maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'DORMQR','LT',m,nrhs,n,-1))

              end if
              if (m >= n) then
                 ! path 1 - overdetermined or exactly determined.
                 maxwrk = max(maxwrk,3*n + (mm + n)*la_ilaenv(1,'DGEBRD',' ',mm,n,-1,- &
                           1))
                 maxwrk = max(maxwrk,3*n + nrhs*la_ilaenv(1,'DORMBR','QLT',mm,nrhs,n,- &
                           1))
                 maxwrk = max(maxwrk,3*n + (n - 1)*la_ilaenv(1,'DORMBR','PLN',n,nrhs,n, &
                           -1))
                 wlalsd = 9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + (smlsiz + 1)**2
                 maxwrk = max(maxwrk,3*n + wlalsd)
                 minwrk = max(3*n + mm,3*n + nrhs,3*n + wlalsd)
              end if
              if (n > m) then
                 wlalsd = 9*m + 2*m*smlsiz + 8*m*nlvl + m*nrhs + (smlsiz + 1)**2
                 if (n >= mnthr) then
                    ! path 2a - underdetermined, with many more columns
                    ! than rows.
                    maxwrk = m + m*la_ilaenv(1,'DGELQF',' ',m,n,-1,-1)
                    maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'DGEBRD',' ',m,m,-1,- &
                              1))
                    maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'DORMBR','QLT',m,nrhs, &
                               m,-1))
                    maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'DORMBR','PLN',m, &
                              nrhs,m,-1))
                    if (nrhs > 1) then
                       maxwrk = max(maxwrk,m*m + m + m*nrhs)
                    else
                       maxwrk = max(maxwrk,m*m + 2*m)
                    end if
                    maxwrk = max(maxwrk,m + nrhs*la_ilaenv(1,'DORMLQ','LT',n,nrhs,m,-1 &
                              ))
                    maxwrk = max(maxwrk,m*m + 4*m + wlalsd)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                    maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                 else
                    ! path 2 - remaining underdetermined cases.
                    maxwrk = 3*m + (n + m)*la_ilaenv(1,'DGEBRD',' ',m,n,-1,-1)
                    maxwrk = max(maxwrk,3*m + nrhs*la_ilaenv(1,'DORMBR','QLT',m,nrhs,n, &
                              -1))
                    maxwrk = max(maxwrk,3*m + m*la_ilaenv(1,'DORMBR','PLN',n,nrhs,m,-1 &
                              ))
                    maxwrk = max(maxwrk,3*m + wlalsd)
                 end if
                 minwrk = max(3*m + nrhs,3*m + m,3*m + wlalsd)
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGELSD',-info)
              return
           else if (lquery) then
              go to 10
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_dlamch('P')
           sfmin = la_dlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_dlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_dlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_dlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_dlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_dlaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_dlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_dlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_dlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure certain entries of b are zero.
           if (m < n) call la_dlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_dgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_dormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_dlaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
                 end if
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              nwork = itaup + n
              ! bidiagonalize r in a.
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_dgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_dormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_dlalsd('U',smlsiz,n,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_dormbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m,wlalsd)) &
                     then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs,4*m + m*lda + &
                        wlalsd)) ldwork = lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_dgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_dlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_dlaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize l in work(il).
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_dgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_dormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_dlalsd('U',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_dormbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_dlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_dormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize a.
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_dgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_dormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_dlalsd('L',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_dormbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_dlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_dlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_dlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_dlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           return
     end subroutine la_dgelsd
#ifdef LA_WITH_XDP
     !> XGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_xgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,iwork, &
               info)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(xdp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,maxmn, &
                     maxwrk,minmn,minwrk,mm,mnthr,nlvl,nwork,smlsiz,wlalsd
           real(xdp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: real,int,log,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           mnthr = la_ilaenv(6,'XGELSD',' ',m,n,nrhs,-1)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           smlsiz = la_ilaenv(9,'XGELSD',' ',0,0,0,0)
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           minwrk = 1
           liwork = 1
           minmn = max(1,minmn)
           nlvl = max(int(log(real(minmn,KIND=xdp)/real(smlsiz + 1,KIND=xdp))/log(two), &
                     KIND=ilp) + 1,0)
           if (info == 0) then
              maxwrk = 0
              liwork = 3*minmn*nlvl + 11*minmn
              mm = m
              if (m >= n .and. m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'XGEQRF',' ',m,n,-1,-1))

                 maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'XORMQR','LT',m,nrhs,n,-1))

              end if
              if (m >= n) then
                 ! path 1 - overdetermined or exactly determined.
                 maxwrk = max(maxwrk,3*n + (mm + n)*la_ilaenv(1,'XGEBRD',' ',mm,n,-1,- &
                           1))
                 maxwrk = max(maxwrk,3*n + nrhs*la_ilaenv(1,'XORMBR','QLT',mm,nrhs,n,- &
                           1))
                 maxwrk = max(maxwrk,3*n + (n - 1)*la_ilaenv(1,'XORMBR','PLN',n,nrhs,n, &
                           -1))
                 wlalsd = 9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + (smlsiz + 1)**2
                 maxwrk = max(maxwrk,3*n + wlalsd)
                 minwrk = max(3*n + mm,3*n + nrhs,3*n + wlalsd)
              end if
              if (n > m) then
                 wlalsd = 9*m + 2*m*smlsiz + 8*m*nlvl + m*nrhs + (smlsiz + 1)**2
                 if (n >= mnthr) then
                    ! path 2a - underdetermined, with many more columns
                    ! than rows.
                    maxwrk = m + m*la_ilaenv(1,'XGELQF',' ',m,n,-1,-1)
                    maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'XGEBRD',' ',m,m,-1,- &
                              1))
                    maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'XORMBR','QLT',m,nrhs, &
                               m,-1))
                    maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'XORMBR','PLN',m, &
                              nrhs,m,-1))
                    if (nrhs > 1) then
                       maxwrk = max(maxwrk,m*m + m + m*nrhs)
                    else
                       maxwrk = max(maxwrk,m*m + 2*m)
                    end if
                    maxwrk = max(maxwrk,m + nrhs*la_ilaenv(1,'XORMLQ','LT',n,nrhs,m,-1 &
                              ))
                    maxwrk = max(maxwrk,m*m + 4*m + wlalsd)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                    maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                 else
                    ! path 2 - remaining underdetermined cases.
                    maxwrk = 3*m + (n + m)*la_ilaenv(1,'XGEBRD',' ',m,n,-1,-1)
                    maxwrk = max(maxwrk,3*m + nrhs*la_ilaenv(1,'XORMBR','QLT',m,nrhs,n, &
                              -1))
                    maxwrk = max(maxwrk,3*m + m*la_ilaenv(1,'XORMBR','PLN',n,nrhs,m,-1 &
                              ))
                    maxwrk = max(maxwrk,3*m + wlalsd)
                 end if
                 minwrk = max(3*m + nrhs,3*m + m,3*m + wlalsd)
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XGELSD',-info)
              return
           else if (lquery) then
              go to 10
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_xlamch('P')
           sfmin = la_xlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_xlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_xlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_xlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_xlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_xlaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_xlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_xlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_xlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure certain entries of b are zero.
           if (m < n) call la_xlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_xgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_xormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_xlaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
                 end if
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              nwork = itaup + n
              ! bidiagonalize r in a.
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_xgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_xormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_xlalsd('U',smlsiz,n,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_xormbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m,wlalsd)) &
                     then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs,4*m + m*lda + &
                        wlalsd)) ldwork = lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_xgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_xlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_xlaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize l in work(il).
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_xgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_xormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_xlalsd('U',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_xormbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_xlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_xormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize a.
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_xgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_xormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_xlalsd('L',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_xormbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_xlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_xlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_xlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_xlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           return
     end subroutine la_xgelsd
#endif
#ifdef LA_WITH_QP
     !> QGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_qgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,iwork, &
               info)
        use la_constants_qp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(qp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,maxmn, &
                     maxwrk,minmn,minwrk,mm,mnthr,nlvl,nwork,smlsiz,wlalsd
           real(qp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: real,int,log,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           mnthr = la_ilaenv(6,'QGELSD',' ',m,n,nrhs,-1)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           smlsiz = la_ilaenv(9,'QGELSD',' ',0,0,0,0)
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           minwrk = 1
           liwork = 1
           minmn = max(1,minmn)
           nlvl = max(int(log(real(minmn,KIND=qp)/real(smlsiz + 1,KIND=qp))/log(two), &
                     KIND=ilp) + 1,0)
           if (info == 0) then
              maxwrk = 0
              liwork = 3*minmn*nlvl + 11*minmn
              mm = m
              if (m >= n .and. m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'QGEQRF',' ',m,n,-1,-1))

                 maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'QORMQR','LT',m,nrhs,n,-1))

              end if
              if (m >= n) then
                 ! path 1 - overdetermined or exactly determined.
                 maxwrk = max(maxwrk,3*n + (mm + n)*la_ilaenv(1,'QGEBRD',' ',mm,n,-1,- &
                           1))
                 maxwrk = max(maxwrk,3*n + nrhs*la_ilaenv(1,'QORMBR','QLT',mm,nrhs,n,- &
                           1))
                 maxwrk = max(maxwrk,3*n + (n - 1)*la_ilaenv(1,'QORMBR','PLN',n,nrhs,n, &
                           -1))
                 wlalsd = 9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + (smlsiz + 1)**2
                 maxwrk = max(maxwrk,3*n + wlalsd)
                 minwrk = max(3*n + mm,3*n + nrhs,3*n + wlalsd)
              end if
              if (n > m) then
                 wlalsd = 9*m + 2*m*smlsiz + 8*m*nlvl + m*nrhs + (smlsiz + 1)**2
                 if (n >= mnthr) then
                    ! path 2a - underdetermined, with many more columns
                    ! than rows.
                    maxwrk = m + m*la_ilaenv(1,'QGELQF',' ',m,n,-1,-1)
                    maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'QGEBRD',' ',m,m,-1,- &
                              1))
                    maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'QORMBR','QLT',m,nrhs, &
                               m,-1))
                    maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'QORMBR','PLN',m, &
                              nrhs,m,-1))
                    if (nrhs > 1) then
                       maxwrk = max(maxwrk,m*m + m + m*nrhs)
                    else
                       maxwrk = max(maxwrk,m*m + 2*m)
                    end if
                    maxwrk = max(maxwrk,m + nrhs*la_ilaenv(1,'QORMLQ','LT',n,nrhs,m,-1 &
                              ))
                    maxwrk = max(maxwrk,m*m + 4*m + wlalsd)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                    maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                 else
                    ! path 2 - remaining underdetermined cases.
                    maxwrk = 3*m + (n + m)*la_ilaenv(1,'QGEBRD',' ',m,n,-1,-1)
                    maxwrk = max(maxwrk,3*m + nrhs*la_ilaenv(1,'QORMBR','QLT',m,nrhs,n, &
                              -1))
                    maxwrk = max(maxwrk,3*m + m*la_ilaenv(1,'QORMBR','PLN',n,nrhs,m,-1 &
                              ))
                    maxwrk = max(maxwrk,3*m + wlalsd)
                 end if
                 minwrk = max(3*m + nrhs,3*m + m,3*m + wlalsd)
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGELSD',-info)
              return
           else if (lquery) then
              go to 10
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_qlamch('P')
           sfmin = la_qlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_qlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_qlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_qlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_qlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_qlaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_qlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_qlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_qlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure certain entries of b are zero.
           if (m < n) call la_qlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns.
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_qgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_qormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_qlaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
                 end if
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              nwork = itaup + n
              ! bidiagonalize r in a.
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_qgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_qormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_qlalsd('U',smlsiz,n,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_qormbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m,wlalsd)) &
                     then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs,4*m + m*lda + &
                        wlalsd)) ldwork = lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_qgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_qlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_qlaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize l in work(il).
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_qgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_qormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_qlalsd('U',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_qormbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_qlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_qormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              nwork = itaup + m
              ! bidiagonalize a.
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_qgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_qormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_qlalsd('L',smlsiz,m,nrhs,s,work(ie),b,ldb,rcond,rank,work( &
                        nwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_qormbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_qlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_qlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_qlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_qlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           return
     end subroutine la_qgelsd
#endif

     !> SGELSS: computes the minimum norm solution to a real linear least
     !> squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_sgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)
        use la_constants_sp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(sp),intent(in) :: rcond
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bdspac,bl,chunk,i,iascl,ibscl,ie,il,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_sgeqrf,lwork_sormqr,lwork_sgebrd,lwork_sormbr,lwork_sorgbr, &
                     lwork_sormlq
           real(sp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'SGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_sgeqrf
                    call la_sgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_sgeqrf = dum(1)
                    ! compute space needed for la_sormqr
                    call la_sormqr('L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_sormqr = dum(1)
                    mm = n
                    maxwrk = max(maxwrk,n + lwork_sgeqrf)
                    maxwrk = max(maxwrk,n + lwork_sormqr)
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute workspace needed for la_sbdsqr
                    bdspac = max(1,5*n)
                    ! compute space needed for la_sgebrd
                    call la_sgebrd(mm,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1,info &
                              )
                    lwork_sgebrd = dum(1)
                    ! compute space needed for la_sormbr
                    call la_sormbr('Q','L','T',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_sormbr = dum(1)
                    ! compute space needed for la_sorgbr
                    call la_sorgbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_sorgbr = dum(1)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,3*n + lwork_sgebrd)
                    maxwrk = max(maxwrk,3*n + lwork_sormbr)
                    maxwrk = max(maxwrk,3*n + lwork_sorgbr)
                    maxwrk = max(maxwrk,bdspac)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = max(3*n + mm,3*n + nrhs,bdspac)
                    maxwrk = max(minwrk,maxwrk)
                 end if
                 if (n > m) then
                    ! compute workspace needed for la_sbdsqr
                    bdspac = max(1,5*m)
                    minwrk = max(3*m + nrhs,3*m + n,bdspac)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_sgebrd
                       call la_sgebrd(m,m,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_sgebrd = dum(1)
                       ! compute space needed for la_sormbr
                       call la_sormbr('Q','L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_sormbr = dum(1)
                       ! compute space needed for la_sorgbr
                       call la_sorgbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_sorgbr = dum(1)
                       ! compute space needed for la_sormlq
                       call la_sormlq('L','T',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_sormlq = dum(1)
                       ! compute total workspace needed
                       maxwrk = m + m*la_ilaenv(1,'SGELQF',' ',m,n,-1,-1)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_sgebrd)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_sormbr)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_sorgbr)
                       maxwrk = max(maxwrk,m*m + m + bdspac)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_sormlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_sgebrd
                       call la_sgebrd(m,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_sgebrd = dum(1)
                       ! compute space needed for la_sormbr
                       call la_sormbr('Q','L','T',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_sormbr = dum(1)
                       ! compute space needed for la_sorgbr
                       call la_sorgbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_sorgbr = dum(1)
                       maxwrk = 3*m + lwork_sgebrd
                       maxwrk = max(maxwrk,3*m + lwork_sormbr)
                       maxwrk = max(maxwrk,3*m + lwork_sorgbr)
                       maxwrk = max(maxwrk,bdspac)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('SGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_slamch('P')
           sfmin = la_slamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_slaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_slaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_slascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_slascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_sgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_sormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_slaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_sgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_sormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (workspace: need 4*n-1, prefer 3*n+(n-1)*nb)
              call la_sorgbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (workspace: need bdspac)
              call la_sbdsqr('U',n,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_srscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_slaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_sgemm('T','N',n,nrhs,n,one,a,lda,b,ldb,zero,work,ldb)

                 call la_slacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_sgemm('T','N',n,bl,n,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_slacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_sgemv('T',n,n,one,a,lda,b,1,zero,work,1)
                 call la_scopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_sgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_slacpy('L',m,m,a,lda,work(il),ldwork)
              call la_slaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_sgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_sormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (workspace: need m*m+5*m-1, prefer m*m+4*m+(m-1)*nb)
              call la_sorgbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of l in work(il) and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need m*m+m+bdspac)
              call la_sbdsqr('U',m,m,0,nrhs,s,work(ie),work(il),ldwork,a,lda,b, &
                         ldb,work(iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_srscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_slaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              iwork = ie
              ! multiply b by right singular vectors of l in work(il)
              ! (workspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_sgemm('T','N',m,nrhs,m,one,work(il),ldwork,b,ldb,zero, &
                           work(iwork),ldb)
                 call la_slacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_sgemm('T','N',m,bl,m,one,work(il),ldwork,b(1,i),ldb, &
                               zero,work(iwork),m)
                    call la_slacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_sgemv('T',m,m,one,work(il),ldwork,b(1,1),1,zero,work( &
                           iwork),1)
                 call la_scopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_slaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_sormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_sgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_sormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (workspace: need 4*m, prefer 3*m+m*nb)
              call la_sorgbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need bdspac)
              call la_sbdsqr('L',m,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_srscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_slaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_sgemm('T','N',n,nrhs,m,one,a,lda,b,ldb,zero,work,ldb)

                 call la_slacpy('F',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_sgemm('T','N',n,bl,m,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_slacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_sgemv('T',m,n,one,a,lda,b,1,zero,work,1)
                 call la_scopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_slascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_slascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_slascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_slascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_sgelss
     !> DGELSS: computes the minimum norm solution to a real linear least
     !> squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_dgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)
        use la_constants_dp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(dp),intent(in) :: rcond
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bdspac,bl,chunk,i,iascl,ibscl,ie,il,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_dgeqrf,lwork_dormqr,lwork_dgebrd,lwork_dormbr,lwork_dorgbr, &
                     lwork_dormlq,lwork_dgelqf
           real(dp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'DGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_dgeqrf
                    call la_dgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_dgeqrf = dum(1)
                    ! compute space needed for la_dormqr
                    call la_dormqr('L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_dormqr = dum(1)
                    mm = n
                    maxwrk = max(maxwrk,n + lwork_dgeqrf)
                    maxwrk = max(maxwrk,n + lwork_dormqr)
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute workspace needed for la_dbdsqr
                    bdspac = max(1,5*n)
                    ! compute space needed for la_dgebrd
                    call la_dgebrd(mm,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1,info &
                              )
                    lwork_dgebrd = dum(1)
                    ! compute space needed for la_dormbr
                    call la_dormbr('Q','L','T',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_dormbr = dum(1)
                    ! compute space needed for la_dorgbr
                    call la_dorgbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_dorgbr = dum(1)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,3*n + lwork_dgebrd)
                    maxwrk = max(maxwrk,3*n + lwork_dormbr)
                    maxwrk = max(maxwrk,3*n + lwork_dorgbr)
                    maxwrk = max(maxwrk,bdspac)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = max(3*n + mm,3*n + nrhs,bdspac)
                    maxwrk = max(minwrk,maxwrk)
                 end if
                 if (n > m) then
                    ! compute workspace needed for la_dbdsqr
                    bdspac = max(1,5*m)
                    minwrk = max(3*m + nrhs,3*m + n,bdspac)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_dgelqf
                       call la_dgelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_dgelqf = dum(1)
                       ! compute space needed for la_dgebrd
                       call la_dgebrd(m,m,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_dgebrd = dum(1)
                       ! compute space needed for la_dormbr
                       call la_dormbr('Q','L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_dormbr = dum(1)
                       ! compute space needed for la_dorgbr
                       call la_dorgbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_dorgbr = dum(1)
                       ! compute space needed for la_dormlq
                       call la_dormlq('L','T',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_dormlq = dum(1)
                       ! compute total workspace needed
                       maxwrk = m + lwork_dgelqf
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_dgebrd)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_dormbr)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_dorgbr)
                       maxwrk = max(maxwrk,m*m + m + bdspac)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_dormlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_dgebrd
                       call la_dgebrd(m,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_dgebrd = dum(1)
                       ! compute space needed for la_dormbr
                       call la_dormbr('Q','L','T',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_dormbr = dum(1)
                       ! compute space needed for la_dorgbr
                       call la_dorgbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_dorgbr = dum(1)
                       maxwrk = 3*m + lwork_dgebrd
                       maxwrk = max(maxwrk,3*m + lwork_dormbr)
                       maxwrk = max(maxwrk,3*m + lwork_dorgbr)
                       maxwrk = max(maxwrk,bdspac)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('DGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_dlamch('P')
           sfmin = la_dlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_dlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_dlaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_dlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_dlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_dgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_dormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_dlaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_dgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_dormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (workspace: need 4*n-1, prefer 3*n+(n-1)*nb)
              call la_dorgbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (workspace: need bdspac)
              call la_dbdsqr('U',n,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_drscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_dlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_dgemm('T','N',n,nrhs,n,one,a,lda,b,ldb,zero,work,ldb)

                 call la_dlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_dgemm('T','N',n,bl,n,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_dlacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_dgemv('T',n,n,one,a,lda,b,1,zero,work,1)
                 call la_dcopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_dgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_dlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_dlaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_dgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_dormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (workspace: need m*m+5*m-1, prefer m*m+4*m+(m-1)*nb)
              call la_dorgbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of l in work(il) and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need m*m+m+bdspac)
              call la_dbdsqr('U',m,m,0,nrhs,s,work(ie),work(il),ldwork,a,lda,b, &
                         ldb,work(iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_drscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_dlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              iwork = ie
              ! multiply b by right singular vectors of l in work(il)
              ! (workspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_dgemm('T','N',m,nrhs,m,one,work(il),ldwork,b,ldb,zero, &
                           work(iwork),ldb)
                 call la_dlacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_dgemm('T','N',m,bl,m,one,work(il),ldwork,b(1,i),ldb, &
                               zero,work(iwork),m)
                    call la_dlacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_dgemv('T',m,m,one,work(il),ldwork,b(1,1),1,zero,work( &
                           iwork),1)
                 call la_dcopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_dlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_dormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_dgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_dormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (workspace: need 4*m, prefer 3*m+m*nb)
              call la_dorgbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need bdspac)
              call la_dbdsqr('L',m,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_drscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_dlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_dgemm('T','N',n,nrhs,m,one,a,lda,b,ldb,zero,work,ldb)

                 call la_dlacpy('F',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_dgemm('T','N',n,bl,m,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_dlacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_dgemv('T',m,n,one,a,lda,b,1,zero,work,1)
                 call la_dcopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_dlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_dlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_dlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_dlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_dgelss
#ifdef LA_WITH_XDP
     !> XGELSS: computes the minimum norm solution to a real linear least
     !> squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_xgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)
        use la_constants_xdp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(xdp),intent(in) :: rcond
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bdspac,bl,chunk,i,iascl,ibscl,ie,il,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_xgeqrf,lwork_xormqr,lwork_xgebrd,lwork_xormbr,lwork_xorgbr, &
                     lwork_xormlq,lwork_dgelqf
           real(xdp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           real(xdp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'XGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_xgeqrf
                    call la_xgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_xgeqrf = dum(1)
                    ! compute space needed for la_xormqr
                    call la_xormqr('L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_xormqr = dum(1)
                    mm = n
                    maxwrk = max(maxwrk,n + lwork_xgeqrf)
                    maxwrk = max(maxwrk,n + lwork_xormqr)
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute workspace needed for la_xbdsqr
                    bdspac = max(1,5*n)
                    ! compute space needed for la_xgebrd
                    call la_xgebrd(mm,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1,info &
                              )
                    lwork_xgebrd = dum(1)
                    ! compute space needed for la_xormbr
                    call la_xormbr('Q','L','T',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_xormbr = dum(1)
                    ! compute space needed for la_xorgbr
                    call la_xorgbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_xorgbr = dum(1)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,3*n + lwork_xgebrd)
                    maxwrk = max(maxwrk,3*n + lwork_xormbr)
                    maxwrk = max(maxwrk,3*n + lwork_xorgbr)
                    maxwrk = max(maxwrk,bdspac)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = max(3*n + mm,3*n + nrhs,bdspac)
                    maxwrk = max(minwrk,maxwrk)
                 end if
                 if (n > m) then
                    ! compute workspace needed for la_xbdsqr
                    bdspac = max(1,5*m)
                    minwrk = max(3*m + nrhs,3*m + n,bdspac)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_xgelqf
                       call la_xgelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_dgelqf = dum(1)
                       ! compute space needed for la_xgebrd
                       call la_xgebrd(m,m,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_xgebrd = dum(1)
                       ! compute space needed for la_xormbr
                       call la_xormbr('Q','L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_xormbr = dum(1)
                       ! compute space needed for la_xorgbr
                       call la_xorgbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_xorgbr = dum(1)
                       ! compute space needed for la_xormlq
                       call la_xormlq('L','T',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_xormlq = dum(1)
                       ! compute total workspace needed
                       maxwrk = m + lwork_dgelqf
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_xgebrd)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_xormbr)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_xorgbr)
                       maxwrk = max(maxwrk,m*m + m + bdspac)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_xormlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_xgebrd
                       call la_xgebrd(m,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_xgebrd = dum(1)
                       ! compute space needed for la_xormbr
                       call la_xormbr('Q','L','T',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_xormbr = dum(1)
                       ! compute space needed for la_xorgbr
                       call la_xorgbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_xorgbr = dum(1)
                       maxwrk = 3*m + lwork_xgebrd
                       maxwrk = max(maxwrk,3*m + lwork_xormbr)
                       maxwrk = max(maxwrk,3*m + lwork_xorgbr)
                       maxwrk = max(maxwrk,bdspac)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('XGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_xlamch('P')
           sfmin = la_xlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_xlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_xlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_xlaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_xlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_xlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_xlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_xgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_xormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_xlaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_xgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_xormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (workspace: need 4*n-1, prefer 3*n+(n-1)*nb)
              call la_xorgbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (workspace: need bdspac)
              call la_xbdsqr('U',n,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_xrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_xlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_xgemm('T','N',n,nrhs,n,one,a,lda,b,ldb,zero,work,ldb)

                 call la_xlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_xgemm('T','N',n,bl,n,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_xlacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_xgemv('T',n,n,one,a,lda,b,1,zero,work,1)
                 call la_xcopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_xgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_xlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_xlaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_xgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_xormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (workspace: need m*m+5*m-1, prefer m*m+4*m+(m-1)*nb)
              call la_xorgbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of l in work(il) and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need m*m+m+bdspac)
              call la_xbdsqr('U',m,m,0,nrhs,s,work(ie),work(il),ldwork,a,lda,b, &
                         ldb,work(iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_xrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_xlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              iwork = ie
              ! multiply b by right singular vectors of l in work(il)
              ! (workspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_xgemm('T','N',m,nrhs,m,one,work(il),ldwork,b,ldb,zero, &
                           work(iwork),ldb)
                 call la_xlacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_xgemm('T','N',m,bl,m,one,work(il),ldwork,b(1,i),ldb, &
                               zero,work(iwork),m)
                    call la_xlacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_xgemv('T',m,m,one,work(il),ldwork,b(1,1),1,zero,work( &
                           iwork),1)
                 call la_xcopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_xlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_xormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_xgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_xormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (workspace: need 4*m, prefer 3*m+m*nb)
              call la_xorgbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need bdspac)
              call la_xbdsqr('L',m,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_xrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_xlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_xgemm('T','N',n,nrhs,m,one,a,lda,b,ldb,zero,work,ldb)

                 call la_xlacpy('F',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_xgemm('T','N',n,bl,m,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_xlacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_xgemv('T',m,n,one,a,lda,b,1,zero,work,1)
                 call la_xcopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_xlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_xlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_xlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_xlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_xgelss
#endif
#ifdef LA_WITH_QP
     !> QGELSS: computes the minimum norm solution to a real linear least
     !> squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_qgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)
        use la_constants_qp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(qp),intent(in) :: rcond
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: s(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bdspac,bl,chunk,i,iascl,ibscl,ie,il,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_qgeqrf,lwork_qormqr,lwork_qgebrd,lwork_qormbr,lwork_qorgbr, &
                     lwork_qormlq,lwork_dgelqf
           real(qp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'QGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_qgeqrf
                    call la_qgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_qgeqrf = dum(1)
                    ! compute space needed for la_qormqr
                    call la_qormqr('L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_qormqr = dum(1)
                    mm = n
                    maxwrk = max(maxwrk,n + lwork_qgeqrf)
                    maxwrk = max(maxwrk,n + lwork_qormqr)
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute workspace needed for la_qbdsqr
                    bdspac = max(1,5*n)
                    ! compute space needed for la_qgebrd
                    call la_qgebrd(mm,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1,info &
                              )
                    lwork_qgebrd = dum(1)
                    ! compute space needed for la_qormbr
                    call la_qormbr('Q','L','T',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_qormbr = dum(1)
                    ! compute space needed for la_qorgbr
                    call la_qorgbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_qorgbr = dum(1)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,3*n + lwork_qgebrd)
                    maxwrk = max(maxwrk,3*n + lwork_qormbr)
                    maxwrk = max(maxwrk,3*n + lwork_qorgbr)
                    maxwrk = max(maxwrk,bdspac)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = max(3*n + mm,3*n + nrhs,bdspac)
                    maxwrk = max(minwrk,maxwrk)
                 end if
                 if (n > m) then
                    ! compute workspace needed for la_qbdsqr
                    bdspac = max(1,5*m)
                    minwrk = max(3*m + nrhs,3*m + n,bdspac)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_qgelqf
                       call la_qgelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_dgelqf = dum(1)
                       ! compute space needed for la_qgebrd
                       call la_qgebrd(m,m,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_qgebrd = dum(1)
                       ! compute space needed for la_qormbr
                       call la_qormbr('Q','L','T',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_qormbr = dum(1)
                       ! compute space needed for la_qorgbr
                       call la_qorgbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_qorgbr = dum(1)
                       ! compute space needed for la_qormlq
                       call la_qormlq('L','T',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_qormlq = dum(1)
                       ! compute total workspace needed
                       maxwrk = m + lwork_dgelqf
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_qgebrd)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_qormbr)
                       maxwrk = max(maxwrk,m*m + 4*m + lwork_qorgbr)
                       maxwrk = max(maxwrk,m*m + m + bdspac)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_qormlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_qgebrd
                       call la_qgebrd(m,n,a,lda,s,dum(1),dum(1),dum(1),dum(1),-1, &
                                 info)
                       lwork_qgebrd = dum(1)
                       ! compute space needed for la_qormbr
                       call la_qormbr('Q','L','T',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_qormbr = dum(1)
                       ! compute space needed for la_qorgbr
                       call la_qorgbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_qorgbr = dum(1)
                       maxwrk = 3*m + lwork_qgebrd
                       maxwrk = max(maxwrk,3*m + lwork_qormbr)
                       maxwrk = max(maxwrk,3*m + lwork_qorgbr)
                       maxwrk = max(maxwrk,bdspac)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('QGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_qlamch('P')
           sfmin = la_qlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',m,n,a,lda,work)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_qlaset('F',max(m,n),nrhs,zero,zero,b,ldb)
              call la_qlaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',m,nrhs,b,ldb,work)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_qlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_qlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (workspace: need 2*n, prefer n+n*nb)
                 call la_qgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (workspace: need n+nrhs, prefer n+nrhs*nb)
                 call la_qormqr('L','T',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_qlaset('L',n - 1,n - 1,zero,zero,a(2,1),lda)
              end if
              ie = 1
              itauq = ie + n
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (workspace: need 3*n+mm, prefer 3*n+(mm+n)*nb)
              call la_qgebrd(mm,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (workspace: need 3*n+nrhs, prefer 3*n+nrhs*nb)
              call la_qormbr('Q','L','T',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (workspace: need 4*n-1, prefer 3*n+(n-1)*nb)
              call la_qorgbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (workspace: need bdspac)
              call la_qbdsqr('U',n,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_qrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_qlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_qgemm('T','N',n,nrhs,n,one,a,lda,b,ldb,zero,work,ldb)

                 call la_qlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_qgemm('T','N',n,bl,n,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_qlacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_qgemv('T',n,n,one,a,lda,b,1,zero,work,1)
                 call la_qcopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (workspace: need 2*m, prefer m+m*nb)
              call la_qgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_qlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_qlaset('U',m - 1,m - 1,zero,zero,work(il + ldwork),ldwork)
              ie = il + ldwork*m
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (workspace: need m*m+5*m, prefer m*m+4*m+2*m*nb)
              call la_qgebrd(m,m,work(il),ldwork,s,work(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (workspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_qormbr('Q','L','T',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (workspace: need m*m+5*m-1, prefer m*m+4*m+(m-1)*nb)
              call la_qorgbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of l in work(il) and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need m*m+m+bdspac)
              call la_qbdsqr('U',m,m,0,nrhs,s,work(ie),work(il),ldwork,a,lda,b, &
                         ldb,work(iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_qrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_qlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              iwork = ie
              ! multiply b by right singular vectors of l in work(il)
              ! (workspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_qgemm('T','N',m,nrhs,m,one,work(il),ldwork,b,ldb,zero, &
                           work(iwork),ldb)
                 call la_qlacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_qgemm('T','N',m,bl,m,one,work(il),ldwork,b(1,i),ldb, &
                               zero,work(iwork),m)
                    call la_qlacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_qgemv('T',m,m,one,work(il),ldwork,b(1,1),1,zero,work( &
                           iwork),1)
                 call la_qcopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_qlaset('F',n - m,nrhs,zero,zero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (workspace: need m+nrhs, prefer m+nrhs*nb)
              call la_qormlq('L','T',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = ie + m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (workspace: need 3*m+n, prefer 3*m+(m+n)*nb)
              call la_qgebrd(m,n,a,lda,s,work(ie),work(itauq),work(itaup),work( &
                        iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (workspace: need 3*m+nrhs, prefer 3*m+nrhs*nb)
              call la_qormbr('Q','L','T',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (workspace: need 4*m, prefer 3*m+m*nb)
              call la_qorgbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              iwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (workspace: need bdspac)
              call la_qbdsqr('L',m,n,0,nrhs,s,work(ie),a,lda,dum,1,b,ldb,work( &
                        iwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_qrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_qlaset('F',1,nrhs,zero,zero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (workspace: need n, prefer n*nrhs)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_qgemm('T','N',n,nrhs,m,one,a,lda,b,ldb,zero,work,ldb)

                 call la_qlacpy('F',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_qgemm('T','N',n,bl,m,one,a,lda,b(1,i),ldb,zero,work, &
                               n)
                    call la_qlacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_qgemv('T',m,n,one,a,lda,b,1,zero,work,1)
                 call la_qcopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_qlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_qlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_qlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_qlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_qgelss
#endif

     !> CGELS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, or its conjugate-transpose, using a QR
     !> or LQ factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**H * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**H * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_cgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_sp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(sp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(sp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'CGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'CUNMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'CUNMQR','LC',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'CGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'CUNMLQ','LC',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'CUNMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=sp)
           end if
           if (info /= 0) then
              call la_xerbla('CGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_claset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_slamch('S')/la_slamch('P')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_clange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_claset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_clange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_cgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
                 call la_cunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,n,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_ctrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**h) * b(1:n,1:nrhs)
                 call la_ctrtrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b, &
                            ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_cunmqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_cgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_ctrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**h * b(1:m,1:nrhs)
                 call la_cunmlq('LEFT','CONJUGATE TRANSPOSE',n,nrhs,m,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**h * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_cunmlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**h) * b(1:m,1:nrhs)
                 call la_ctrtrs('LOWER','CONJUGATE TRANSPOSE','NON-UNIT',m,nrhs,a,lda, &
                           b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_clascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_clascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_clascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_clascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=sp)
           return
     end subroutine la_cgels
     !> ZGELS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, or its conjugate-transpose, using a QR
     !> or LQ factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**H * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**H * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_zgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_dp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(dp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(dp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'ZGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'ZUNMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'ZUNMQR','LC',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'ZGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'ZUNMLQ','LC',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'ZUNMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=dp)
           end if
           if (info /= 0) then
              call la_xerbla('ZGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_zlaset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_dlamch('S')/la_dlamch('P')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_zlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_zlange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_zgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
                 call la_zunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,n,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_ztrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**h) * b(1:n,1:nrhs)
                 call la_ztrtrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b, &
                            ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_zunmqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_zgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_ztrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**h * b(1:m,1:nrhs)
                 call la_zunmlq('LEFT','CONJUGATE TRANSPOSE',n,nrhs,m,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**h * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_zunmlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**h) * b(1:m,1:nrhs)
                 call la_ztrtrs('LOWER','CONJUGATE TRANSPOSE','NON-UNIT',m,nrhs,a,lda, &
                           b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_zlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_zlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_zlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_zlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=dp)
           return
     end subroutine la_zgels
#ifdef LA_WITH_XDP
     !> YGELS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, or its conjugate-transpose, using a QR
     !> or LQ factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**H * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**H * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_ygels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_xdp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(xdp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(xdp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'YGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'YUNMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'YUNMQR','LC',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'YGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'YUNMLQ','LC',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'YUNMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=xdp)
           end if
           if (info /= 0) then
              call la_xerbla('YGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_ylaset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_xlamch('S')/la_xlamch('P')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_ylange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_ylaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_ylange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_ygeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
                 call la_yunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,n,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_ytrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**h) * b(1:n,1:nrhs)
                 call la_ytrtrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b, &
                            ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_yunmqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_ygelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_ytrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**h * b(1:m,1:nrhs)
                 call la_yunmlq('LEFT','CONJUGATE TRANSPOSE',n,nrhs,m,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**h * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_yunmlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**h) * b(1:m,1:nrhs)
                 call la_ytrtrs('LOWER','CONJUGATE TRANSPOSE','NON-UNIT',m,nrhs,a,lda, &
                           b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_ylascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_ylascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_ylascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_ylascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=xdp)
           return
     end subroutine la_ygels
#endif
#ifdef LA_WITH_QP
     !> WGELS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, or its conjugate-transpose, using a QR
     !> or LQ factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an underdetermined system A**H * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**H * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_wgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_qp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tpsd
           integer(ilp) :: brow,i,iascl,ibscl,j,mn,nb,scllen,wsize
           real(qp) :: anrm,bignum,bnrm,smlnum
           ! Local Arrays
           real(qp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input arguments.
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           else if (lwork < max(1,mn + max(mn,nrhs)) .and. .not. lquery) then
              info = -10
           end if
           ! figure out optimal block size
           if (info == 0 .or. info == -10) then
              tpsd = .true.
              if (la_lsame(trans,'N')) tpsd = .false.
              if (m >= n) then
                 nb = la_ilaenv(1,'WGEQRF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'WUNMQR','LN',m,nrhs,n,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'WUNMQR','LC',m,nrhs,n,-1))
                 end if
              else
                 nb = la_ilaenv(1,'WGELQF',' ',m,n,-1,-1)
                 if (tpsd) then
                    nb = max(nb,la_ilaenv(1,'WUNMLQ','LC',n,nrhs,m,-1))
                 else
                    nb = max(nb,la_ilaenv(1,'WUNMLQ','LN',n,nrhs,m,-1))
                 end if
              end if
              wsize = max(1,mn + max(mn,nrhs)*nb)
              work(1) = real(wsize,KIND=qp)
           end if
           if (info /= 0) then
              call la_xerbla('WGELS ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              call la_wlaset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
              return
           end if
           ! get machine parameters
           smlnum = la_qlamch('S')/la_qlamch('P')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_wlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tpsd) brow = n
           bnrm = la_wlange('M',brow,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
              call la_wgeqrf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least n, optimally n*nb
              if (.not. tpsd) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
                 call la_wunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,n,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
                 call la_wtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 scllen = n
              else
                 ! underdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**h) * b(1:n,1:nrhs)
                 call la_wtrtrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT',n,nrhs,a,lda,b, &
                            ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = zero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_wunmqr('LEFT','NO TRANSPOSE',m,nrhs,n,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_wgelqf(m,n,a,lda,work(1),work(mn + 1),lwork - mn,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tpsd) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_wtrtrs('LOWER','NO TRANSPOSE','NON-UNIT',m,nrhs,a,lda,b,ldb, &
                           info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**h * b(1:m,1:nrhs)
                 call la_wunmlq('LEFT','CONJUGATE TRANSPOSE',n,nrhs,m,a,lda,work(1), &
                           b,ldb,work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**h * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_wunmlq('LEFT','NO TRANSPOSE',n,nrhs,m,a,lda,work(1),b,ldb, &
                            work(mn + 1),lwork - mn,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**h) * b(1:m,1:nrhs)
                 call la_wtrtrs('LOWER','CONJUGATE TRANSPOSE','NON-UNIT',m,nrhs,a,lda, &
                           b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_wlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
              call la_wlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
              call la_wlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_wlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(wsize,KIND=qp)
           return
     end subroutine la_wgels
#endif

     !> CGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_cgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               iwork,info)
        use la_constants_sp,only:zero,one,two,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(sp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: rwork(*),s(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,lrwork, &
                     maxmn,maxwrk,minmn,minwrk,mm,mnthr,nlvl,nrwork,nwork,smlsiz
           real(sp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: int,log,max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              liwork = 1
              lrwork = 1
              if (minmn > 0) then
                 smlsiz = la_ilaenv(9,'CGELSD',' ',0,0,0,0)
                 mnthr = la_ilaenv(6,'CGELSD',' ',m,n,nrhs,-1)
                 nlvl = max(int(log(real(minmn,KIND=sp)/real(smlsiz + 1,KIND=sp))/log( &
                           two),KIND=ilp) + 1,0)
                 liwork = 3*minmn*nlvl + 11*minmn
                 mm = m
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns.
                    mm = n
                    maxwrk = max(maxwrk,n*la_ilaenv(1,'CGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,nrhs*la_ilaenv(1,'CUNMQR','LC',m,nrhs,n,-1))

                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined.
                    lrwork = 10*n + 2*n*smlsiz + 8*n*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    maxwrk = max(maxwrk,2*n + (mm + n)*la_ilaenv(1,'CGEBRD',' ',mm,n, &
                              -1,-1))
                    maxwrk = max(maxwrk,2*n + nrhs*la_ilaenv(1,'CUNMBR','QLC',mm,nrhs, &
                              n,-1))
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'CUNMBR','PLN',n, &
                              nrhs,n,-1))
                    maxwrk = max(maxwrk,2*n + n*nrhs)
                    minwrk = max(2*n + mm,2*n + n*nrhs)
                 end if
                 if (n > m) then
                    lrwork = 10*m + 2*m*smlsiz + 8*m*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                                 ! than rows.
                       maxwrk = m + m*la_ilaenv(1,'CGELQF',' ',m,n,-1,-1)
                       maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'CGEBRD',' ',m,m, &
                                  -1,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'CUNMBR','QLC',m, &
                                  nrhs,m,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'CUNMLQ', &
                                 'LC',n,nrhs,m,-1))
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m*m + 4*m + m*nrhs)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                       maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                    else
                       ! path 2 - underdetermined.
                       maxwrk = 2*m + (n + m)*la_ilaenv(1,'CGEBRD',' ',m,n,-1,-1)

                       maxwrk = max(maxwrk,2*m + nrhs*la_ilaenv(1,'CUNMBR','QLC',m,nrhs, &
                                  m,-1))
                       maxwrk = max(maxwrk,2*m + m*la_ilaenv(1,'CUNMBR','PLN',n,nrhs,m, &
                                  -1))
                       maxwrk = max(maxwrk,2*m + m*nrhs)
                    end if
                    minwrk = max(2*m + n,2*m + m*nrhs)
                 end if
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              rwork(1) = lrwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGELSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_slamch('P')
           sfmin = la_slamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_clange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_clascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_claset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_slaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_clange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_clascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_clascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure b(m+1:n,:) = 0
           if (m < n) call la_claset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (rworkspace: need n)
                 ! (cworkspace: need n, prefer n*nb)
                 call la_cgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (rworkspace: need n)
                 ! (cworkspace: need nrhs, prefer nrhs*nb)
                 call la_cunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_claset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
                 end if
              end if
              itauq = 1
              itaup = itauq + n
              nwork = itaup + n
              ie = 1
              nrwork = ie + n
              ! bidiagonalize r in a.
              ! (rworkspace: need n)
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              call la_cgebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              call la_cunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_clalsd('U',smlsiz,n,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_cunmbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (cworkspace: need 2*m, prefer m+m*nb)
              call la_cgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_clacpy('L',m,m,a,lda,work(il),ldwork)
              call la_claset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              itauq = il + ldwork*m
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize l in work(il).
              ! (rworkspace: need m)
              ! (cworkspace: need m*m+4*m, prefer m*m+4*m+2*m*nb)
              call la_cgebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (cworkspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_cunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_clalsd('U',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_cunmbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_claset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (cworkspace: need nrhs, prefer nrhs*nb)
              call la_cunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              itauq = 1
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize a.
              ! (rworkspace: need m)
              ! (cworkspace: need 2*m+n, prefer 2*m+(m+n)*nb)
              call la_cgebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              call la_cunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_clalsd('L',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_cunmbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_clascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_clascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_clascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_clascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           rwork(1) = lrwork
           return
     end subroutine la_cgelsd
     !> ZGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_zgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               iwork,info)
        use la_constants_dp,only:zero,one,two,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(dp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: rwork(*),s(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,lrwork, &
                     maxmn,maxwrk,minmn,minwrk,mm,mnthr,nlvl,nrwork,nwork,smlsiz
           real(dp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: int,log,max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              liwork = 1
              lrwork = 1
              if (minmn > 0) then
                 smlsiz = la_ilaenv(9,'ZGELSD',' ',0,0,0,0)
                 mnthr = la_ilaenv(6,'ZGELSD',' ',m,n,nrhs,-1)
                 nlvl = max(int(log(real(minmn,KIND=dp)/real(smlsiz + 1,KIND=dp))/log( &
                           two),KIND=ilp) + 1,0)
                 liwork = 3*minmn*nlvl + 11*minmn
                 mm = m
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns.
                    mm = n
                    maxwrk = max(maxwrk,n*la_ilaenv(1,'ZGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,nrhs*la_ilaenv(1,'ZUNMQR','LC',m,nrhs,n,-1))

                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined.
                    lrwork = 10*n + 2*n*smlsiz + 8*n*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    maxwrk = max(maxwrk,2*n + (mm + n)*la_ilaenv(1,'ZGEBRD',' ',mm,n, &
                              -1,-1))
                    maxwrk = max(maxwrk,2*n + nrhs*la_ilaenv(1,'ZUNMBR','QLC',mm,nrhs, &
                              n,-1))
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'ZUNMBR','PLN',n, &
                              nrhs,n,-1))
                    maxwrk = max(maxwrk,2*n + n*nrhs)
                    minwrk = max(2*n + mm,2*n + n*nrhs)
                 end if
                 if (n > m) then
                    lrwork = 10*m + 2*m*smlsiz + 8*m*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                                 ! than rows.
                       maxwrk = m + m*la_ilaenv(1,'ZGELQF',' ',m,n,-1,-1)
                       maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'ZGEBRD',' ',m,m, &
                                  -1,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'ZUNMBR','QLC',m, &
                                  nrhs,m,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'ZUNMLQ', &
                                 'LC',n,nrhs,m,-1))
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m*m + 4*m + m*nrhs)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                       maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                    else
                       ! path 2 - underdetermined.
                       maxwrk = 2*m + (n + m)*la_ilaenv(1,'ZGEBRD',' ',m,n,-1,-1)

                       maxwrk = max(maxwrk,2*m + nrhs*la_ilaenv(1,'ZUNMBR','QLC',m,nrhs, &
                                  m,-1))
                       maxwrk = max(maxwrk,2*m + m*la_ilaenv(1,'ZUNMBR','PLN',n,nrhs,m, &
                                  -1))
                       maxwrk = max(maxwrk,2*m + m*nrhs)
                    end if
                    minwrk = max(2*m + n,2*m + m*nrhs)
                 end if
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              rwork(1) = lrwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGELSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_dlamch('P')
           sfmin = la_dlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_zlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_zlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_zlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_dlaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_zlange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_zlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_zlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure b(m+1:n,:) = 0
           if (m < n) call la_zlaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (rworkspace: need n)
                 ! (cworkspace: need n, prefer n*nb)
                 call la_zgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (rworkspace: need n)
                 ! (cworkspace: need nrhs, prefer nrhs*nb)
                 call la_zunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_zlaset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
                 end if
              end if
              itauq = 1
              itaup = itauq + n
              nwork = itaup + n
              ie = 1
              nrwork = ie + n
              ! bidiagonalize r in a.
              ! (rworkspace: need n)
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              call la_zgebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              call la_zunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_zlalsd('U',smlsiz,n,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_zunmbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (cworkspace: need 2*m, prefer m+m*nb)
              call la_zgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_zlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_zlaset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              itauq = il + ldwork*m
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize l in work(il).
              ! (rworkspace: need m)
              ! (cworkspace: need m*m+4*m, prefer m*m+4*m+2*m*nb)
              call la_zgebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (cworkspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_zunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_zlalsd('U',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_zunmbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_zlaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (cworkspace: need nrhs, prefer nrhs*nb)
              call la_zunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              itauq = 1
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize a.
              ! (rworkspace: need m)
              ! (cworkspace: need 2*m+n, prefer 2*m+(m+n)*nb)
              call la_zgebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              call la_zunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_zlalsd('L',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_zunmbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_zlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_zlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_zlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_zlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           rwork(1) = lrwork
           return
     end subroutine la_zgelsd
#ifdef LA_WITH_XDP
     !> YGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_ygelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               iwork,info)
        use la_constants_xdp,only:zero,one,two,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(xdp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(out) :: rwork(*),s(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,lrwork, &
                     maxmn,maxwrk,minmn,minwrk,mm,mnthr,nlvl,nrwork,nwork,smlsiz
           real(xdp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: int,log,max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              liwork = 1
              lrwork = 1
              if (minmn > 0) then
                 smlsiz = la_ilaenv(9,'YGELSD',' ',0,0,0,0)
                 mnthr = la_ilaenv(6,'YGELSD',' ',m,n,nrhs,-1)
                 nlvl = max(int(log(real(minmn,KIND=xdp)/real(smlsiz + 1,KIND=xdp))/log( &
                           two),KIND=ilp) + 1,0)
                 liwork = 3*minmn*nlvl + 11*minmn
                 mm = m
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns.
                    mm = n
                    maxwrk = max(maxwrk,n*la_ilaenv(1,'YGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,nrhs*la_ilaenv(1,'YUNMQR','LC',m,nrhs,n,-1))

                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined.
                    lrwork = 10*n + 2*n*smlsiz + 8*n*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    maxwrk = max(maxwrk,2*n + (mm + n)*la_ilaenv(1,'YGEBRD',' ',mm,n, &
                              -1,-1))
                    maxwrk = max(maxwrk,2*n + nrhs*la_ilaenv(1,'YUNMBR','QLC',mm,nrhs, &
                              n,-1))
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'YUNMBR','PLN',n, &
                              nrhs,n,-1))
                    maxwrk = max(maxwrk,2*n + n*nrhs)
                    minwrk = max(2*n + mm,2*n + n*nrhs)
                 end if
                 if (n > m) then
                    lrwork = 10*m + 2*m*smlsiz + 8*m*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                                 ! than rows.
                       maxwrk = m + m*la_ilaenv(1,'YGELQF',' ',m,n,-1,-1)
                       maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'YGEBRD',' ',m,m, &
                                  -1,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'YUNMBR','QLC',m, &
                                  nrhs,m,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'YUNMLQ', &
                                 'LC',n,nrhs,m,-1))
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m*m + 4*m + m*nrhs)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                       maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                    else
                       ! path 2 - underdetermined.
                       maxwrk = 2*m + (n + m)*la_ilaenv(1,'YGEBRD',' ',m,n,-1,-1)

                       maxwrk = max(maxwrk,2*m + nrhs*la_ilaenv(1,'YUNMBR','QLC',m,nrhs, &
                                  m,-1))
                       maxwrk = max(maxwrk,2*m + m*la_ilaenv(1,'YUNMBR','PLN',n,nrhs,m, &
                                  -1))
                       maxwrk = max(maxwrk,2*m + m*nrhs)
                    end if
                    minwrk = max(2*m + n,2*m + m*nrhs)
                 end if
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              rwork(1) = lrwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YGELSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_xlamch('P')
           sfmin = la_xlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_ylange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_ylascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_ylaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_xlaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_ylange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_ylascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_ylascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure b(m+1:n,:) = 0
           if (m < n) call la_ylaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (rworkspace: need n)
                 ! (cworkspace: need n, prefer n*nb)
                 call la_ygeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (rworkspace: need n)
                 ! (cworkspace: need nrhs, prefer nrhs*nb)
                 call la_yunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_ylaset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
                 end if
              end if
              itauq = 1
              itaup = itauq + n
              nwork = itaup + n
              ie = 1
              nrwork = ie + n
              ! bidiagonalize r in a.
              ! (rworkspace: need n)
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              call la_ygebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              call la_yunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_ylalsd('U',smlsiz,n,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_yunmbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (cworkspace: need 2*m, prefer m+m*nb)
              call la_ygelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_ylacpy('L',m,m,a,lda,work(il),ldwork)
              call la_ylaset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              itauq = il + ldwork*m
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize l in work(il).
              ! (rworkspace: need m)
              ! (cworkspace: need m*m+4*m, prefer m*m+4*m+2*m*nb)
              call la_ygebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (cworkspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_yunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_ylalsd('U',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_yunmbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_ylaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (cworkspace: need nrhs, prefer nrhs*nb)
              call la_yunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              itauq = 1
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize a.
              ! (rworkspace: need m)
              ! (cworkspace: need 2*m+n, prefer 2*m+(m+n)*nb)
              call la_ygebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              call la_yunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_ylalsd('L',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_yunmbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_ylascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_ylascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_ylascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_ylascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           rwork(1) = lrwork
           return
     end subroutine la_ygelsd
#endif
#ifdef LA_WITH_QP
     !> WGELSD: computes the minimum-norm solution to a real linear least
     !> squares problem:
     !> minimize 2-norm(| b - A*x |)
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The problem is solved in three steps:
     !> (1) Reduce the coefficient matrix A to bidiagonal form with
     !> Householder transformations, reducing the original problem
     !> into a "bidiagonal least squares problem" (BLS)
     !> (2) Solve the BLS using a divide and conquer approach.
     !> (3) Apply back all the Householder transformations to solve
     !> the original least squares problem.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     subroutine la_wgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               iwork,info)
        use la_constants_qp,only:zero,one,two,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(qp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: rwork(*),s(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: iascl,ibscl,ie,il,itau,itaup,itauq,ldwork,liwork,lrwork, &
                     maxmn,maxwrk,minmn,minwrk,mm,mnthr,nlvl,nrwork,nwork,smlsiz
           real(qp) :: anrm,bignum,bnrm,eps,sfmin,smlnum
           ! Intrinsic Functions
           intrinsic :: int,log,max,min,real
           ! Executable Statements
           ! test the input arguments.
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace.
           ! (note: comments in the code beginning "workspace:" describe the
           ! minimal amount of workspace needed at that point in the code,
           ! as well as the preferred amount for good performance.
           ! nb refers to the optimal block size for the immediately
           ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              liwork = 1
              lrwork = 1
              if (minmn > 0) then
                 smlsiz = la_ilaenv(9,'WGELSD',' ',0,0,0,0)
                 mnthr = la_ilaenv(6,'WGELSD',' ',m,n,nrhs,-1)
                 nlvl = max(int(log(real(minmn,KIND=qp)/real(smlsiz + 1,KIND=qp))/log( &
                           two),KIND=ilp) + 1,0)
                 liwork = 3*minmn*nlvl + 11*minmn
                 mm = m
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns.
                    mm = n
                    maxwrk = max(maxwrk,n*la_ilaenv(1,'WGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,nrhs*la_ilaenv(1,'WUNMQR','LC',m,nrhs,n,-1))

                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined.
                    lrwork = 10*n + 2*n*smlsiz + 8*n*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    maxwrk = max(maxwrk,2*n + (mm + n)*la_ilaenv(1,'WGEBRD',' ',mm,n, &
                              -1,-1))
                    maxwrk = max(maxwrk,2*n + nrhs*la_ilaenv(1,'WUNMBR','QLC',mm,nrhs, &
                              n,-1))
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'WUNMBR','PLN',n, &
                              nrhs,n,-1))
                    maxwrk = max(maxwrk,2*n + n*nrhs)
                    minwrk = max(2*n + mm,2*n + n*nrhs)
                 end if
                 if (n > m) then
                    lrwork = 10*m + 2*m*smlsiz + 8*m*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*( &
                              1 + nrhs) + 2*nrhs)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                                 ! than rows.
                       maxwrk = m + m*la_ilaenv(1,'WGELQF',' ',m,n,-1,-1)
                       maxwrk = max(maxwrk,m*m + 4*m + 2*m*la_ilaenv(1,'WGEBRD',' ',m,m, &
                                  -1,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + nrhs*la_ilaenv(1,'WUNMBR','QLC',m, &
                                  nrhs,m,-1))
                       maxwrk = max(maxwrk,m*m + 4*m + (m - 1)*la_ilaenv(1,'WUNMLQ', &
                                 'LC',n,nrhs,m,-1))
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m*m + 4*m + m*nrhs)
           ! xxx: ensure the path 2a case below is triggered.  the workspace
           ! calculation should use queries for all routines eventually.
                       maxwrk = max(maxwrk,4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m))
                    else
                       ! path 2 - underdetermined.
                       maxwrk = 2*m + (n + m)*la_ilaenv(1,'WGEBRD',' ',m,n,-1,-1)

                       maxwrk = max(maxwrk,2*m + nrhs*la_ilaenv(1,'WUNMBR','QLC',m,nrhs, &
                                  m,-1))
                       maxwrk = max(maxwrk,2*m + m*la_ilaenv(1,'WUNMBR','PLN',n,nrhs,m, &
                                  -1))
                       maxwrk = max(maxwrk,2*m + m*nrhs)
                    end if
                    minwrk = max(2*m + n,2*m + m*nrhs)
                 end if
              end if
              minwrk = min(minwrk,maxwrk)
              work(1) = maxwrk
              iwork(1) = liwork
              rwork(1) = lrwork
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGELSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters.
           eps = la_qlamch('P')
           sfmin = la_qlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a if max entry outside range [smlnum,bignum].
           anrm = la_wlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_wlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_wlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_qlaset('F',minmn,1,zero,zero,s,1)
              rank = 0
              go to 10
           end if
           ! scale b if max entry outside range [smlnum,bignum].
           bnrm = la_wlange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum.
              call la_wlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum.
              call la_wlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! if m < n make sure b(m+1:n,:) = 0
           if (m < n) call la_wlaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
           ! overdetermined case.
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined.
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 nwork = itau + n
                 ! compute a=q*r.
                 ! (rworkspace: need n)
                 ! (cworkspace: need n, prefer n*nb)
                 call la_wgeqrf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1, &
                           info)
                 ! multiply b by transpose(q).
                 ! (rworkspace: need n)
                 ! (cworkspace: need nrhs, prefer nrhs*nb)
                 call la_wunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           nwork),lwork - nwork + 1,info)
                 ! zero out below r.
                 if (n > 1) then
                    call la_wlaset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
                 end if
              end if
              itauq = 1
              itaup = itauq + n
              nwork = itaup + n
              ie = 1
              nrwork = ie + n
              ! bidiagonalize r in a.
              ! (rworkspace: need n)
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              call la_wgebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r.
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              call la_wunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_wlalsd('U',smlsiz,n,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of r.
              call la_wunmbr('P','L','N',n,nrhs,n,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           else if (n >= mnthr .and. lwork >= 4*m + m*m + max(m,2*m - 4,nrhs,n - 3*m)) then
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm.
              ldwork = m
              if (lwork >= max(4*m + m*lda + max(m,2*m - 4,nrhs,n - 3*m),m*lda + m + m*nrhs)) ldwork = &
                        lda
              itau = 1
              nwork = m + 1
              ! compute a=l*q.
              ! (cworkspace: need 2*m, prefer m+m*nb)
              call la_wgelqf(m,n,a,lda,work(itau),work(nwork),lwork - nwork + 1,info)

              il = nwork
              ! copy l to work(il), zeroing out above its diagonal.
              call la_wlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_wlaset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              itauq = il + ldwork*m
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize l in work(il).
              ! (rworkspace: need m)
              ! (cworkspace: need m*m+4*m, prefer m*m+4*m+2*m*nb)
              call la_wgebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l.
              ! (cworkspace: need m*m+4*m+nrhs, prefer m*m+4*m+nrhs*nb)
              call la_wunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_wlalsd('U',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of l.
              call la_wunmbr('P','L','N',m,nrhs,m,work(il),ldwork,work(itaup),b, &
                        ldb,work(nwork),lwork - nwork + 1,info)
              ! zero out below first m rows of b.
              call la_wlaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              nwork = itau + m
              ! multiply transpose(q) by b.
              ! (cworkspace: need nrhs, prefer nrhs*nb)
              call la_wunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(nwork) &
                        ,lwork - nwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases.
              itauq = 1
              itaup = itauq + m
              nwork = itaup + m
              ie = 1
              nrwork = ie + m
              ! bidiagonalize a.
              ! (rworkspace: need m)
              ! (cworkspace: need 2*m+n, prefer 2*m+(m+n)*nb)
              call la_wgebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         nwork),lwork - nwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors.
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              call la_wunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
              ! solve the bidiagonal least squares problem.
              call la_wlalsd('L',smlsiz,m,nrhs,s,rwork(ie),b,ldb,rcond,rank,work( &
                        nwork),rwork(nrwork),iwork,info)
              if (info /= 0) then
                 go to 10
              end if
              ! multiply b by right bidiagonalizing vectors of a.
              call la_wunmbr('P','L','N',n,nrhs,m,a,lda,work(itaup),b,ldb,work( &
                        nwork),lwork - nwork + 1,info)
           end if
           ! undo scaling.
           if (iascl == 1) then
              call la_wlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_wlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_wlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_wlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           10 continue
           work(1) = maxwrk
           iwork(1) = liwork
           rwork(1) = lrwork
           return
     end subroutine la_wgelsd
#endif

     !> CGELSS: computes the minimum norm solution to a complex linear
     !> least squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_cgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(sp),intent(in) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: rwork(*),s(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bl,chunk,i,iascl,ibscl,ie,il,irwork,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_cgeqrf,lwork_cunmqr,lwork_cgebrd,lwork_cunmbr,lwork_cungbr, &
                     lwork_cunmlq,lwork_cgelqf
           real(sp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           complex(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace refers
             ! to real workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'CGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_cgeqrf
                    call la_cgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_cgeqrf = real(dum(1),KIND=sp)
                    ! compute space needed for la_cunmqr
                    call la_cunmqr('L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_cunmqr = real(dum(1),KIND=sp)
                    mm = n
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'CGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'CUNMQR','LC',m,nrhs,n,- &
                              1))
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute space needed for la_cgebrd
                    call la_cgebrd(mm,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                    lwork_cgebrd = real(dum(1),KIND=sp)
                    ! compute space needed for la_cunmbr
                    call la_cunmbr('Q','L','C',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_cunmbr = real(dum(1),KIND=sp)
                    ! compute space needed for la_cungbr
                    call la_cungbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_cungbr = real(dum(1),KIND=sp)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,2*n + lwork_cgebrd)
                    maxwrk = max(maxwrk,2*n + lwork_cunmbr)
                    maxwrk = max(maxwrk,2*n + lwork_cungbr)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = 2*n + max(nrhs,m)
                 end if
                 if (n > m) then
                    minwrk = 2*m + max(nrhs,n)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_cgelqf
                       call la_cgelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_cgelqf = real(dum(1),KIND=sp)
                       ! compute space needed for la_cgebrd
                       call la_cgebrd(m,m,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_cgebrd = real(dum(1),KIND=sp)
                       ! compute space needed for la_cunmbr
                       call la_cunmbr('Q','L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_cunmbr = real(dum(1),KIND=sp)
                       ! compute space needed for la_cungbr
                       call la_cungbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_cungbr = real(dum(1),KIND=sp)
                       ! compute space needed for la_cunmlq
                       call la_cunmlq('L','C',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_cunmlq = real(dum(1),KIND=sp)
                       ! compute total workspace needed
                       maxwrk = m + lwork_cgelqf
                       maxwrk = max(maxwrk,3*m + m*m + lwork_cgebrd)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_cunmbr)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_cungbr)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_cunmlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_cgebrd
                       call la_cgebrd(m,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_cgebrd = real(dum(1),KIND=sp)
                       ! compute space needed for la_cunmbr
                       call la_cunmbr('Q','L','C',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_cunmbr = real(dum(1),KIND=sp)
                       ! compute space needed for la_cungbr
                       call la_cungbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_cungbr = real(dum(1),KIND=sp)
                       maxwrk = 2*m + lwork_cgebrd
                       maxwrk = max(maxwrk,2*m + lwork_cunmbr)
                       maxwrk = max(maxwrk,2*m + lwork_cungbr)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('CGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_slamch('P')
           sfmin = la_slamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_claset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_slaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (cworkspace: need 2*n, prefer n+n*nb)
                 ! (rworkspace: none)
                 call la_cgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (cworkspace: need n+nrhs, prefer n+nrhs*nb)
                 ! (rworkspace: none)
                 call la_cunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_claset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
              end if
              ie = 1
              itauq = 1
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              ! (rworkspace: need n)
              call la_cgebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              ! (rworkspace: none)
              call la_cunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (cworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cungbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_cbdsqr('U',n,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_csrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_claset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_cgemm('C','N',n,nrhs,n,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_clacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_cgemm('C','N',n,bl,n,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_clacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_cgemv('C',n,n,cone,a,lda,b,1,czero,work,1)
                 call la_ccopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 3*m + m*m + max(m,nrhs,n - 2*m)) then
              ! underdetermined case, m much less than n
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= 3*m + m*lda + max(m,nrhs,n - 2*m)) ldwork = lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (cworkspace: need 2*m, prefer m+m*nb)
              ! (rworkspace: none)
              call la_cgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_clacpy('L',m,m,a,lda,work(il),ldwork)
              call la_claset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              ie = 1
              itauq = il + ldwork*m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (cworkspace: need m*m+4*m, prefer m*m+3*m+2*m*nb)
              ! (rworkspace: need m)
              call la_cgebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (cworkspace: need m*m+3*m+nrhs, prefer m*m+3*m+nrhs*nb)
              ! (rworkspace: none)
              call la_cunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (cworkspace: need m*m+4*m-1, prefer m*m+3*m+(m-1)*nb)
              ! (rworkspace: none)
              call la_cungbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration, computing right singular
              ! vectors of l in work(il) and multiplying b by transpose of
              ! left singular vectors
              ! (cworkspace: need m*m)
              ! (rworkspace: need bdspac)
              call la_cbdsqr('U',m,m,0,nrhs,s,rwork(ie),work(il),ldwork,a,lda, &
                        b,ldb,rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_csrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_claset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              iwork = il + m*ldwork
              ! multiply b by right singular vectors of l in work(il)
              ! (cworkspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_cgemm('C','N',m,nrhs,m,cone,work(il),ldwork,b,ldb,czero, &
                           work(iwork),ldb)
                 call la_clacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_cgemm('C','N',m,bl,m,cone,work(il),ldwork,b(1,i), &
                              ldb,czero,work(iwork),m)
                    call la_clacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_cgemv('C',m,m,cone,work(il),ldwork,b(1,1),1,czero,work( &
                            iwork),1)
                 call la_ccopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_claset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (cworkspace: need m+nrhs, prefer m+nhrs*nb)
              ! (rworkspace: none)
              call la_cunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = 1
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (cworkspace: need 3*m, prefer 2*m+(m+n)*nb)
              ! (rworkspace: need n)
              call la_cgebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              ! (rworkspace: none)
              call la_cunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (cworkspace: need 3*m, prefer 2*m+m*nb)
              ! (rworkspace: none)
              call la_cungbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_cbdsqr('L',m,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_csrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_claset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_cgemm('C','N',n,nrhs,m,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_clacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_cgemm('C','N',n,bl,m,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_clacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_cgemv('C',m,n,cone,a,lda,b,1,czero,work,1)
                 call la_ccopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_clascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_clascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_slascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_clascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_clascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_cgelss
     !> ZGELSS: computes the minimum norm solution to a complex linear
     !> least squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_zgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(dp),intent(in) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: rwork(*),s(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bl,chunk,i,iascl,ibscl,ie,il,irwork,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_zgeqrf,lwork_zunmqr,lwork_zgebrd,lwork_zunmbr,lwork_zungbr, &
                     lwork_zunmlq,lwork_zgelqf
           real(dp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           complex(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace refers
             ! to real workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'ZGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_zgeqrf
                    call la_zgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_zgeqrf = real(dum(1),KIND=dp)
                    ! compute space needed for la_zunmqr
                    call la_zunmqr('L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_zunmqr = real(dum(1),KIND=dp)
                    mm = n
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'ZGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'ZUNMQR','LC',m,nrhs,n,- &
                              1))
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute space needed for la_zgebrd
                    call la_zgebrd(mm,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                    lwork_zgebrd = real(dum(1),KIND=dp)
                    ! compute space needed for la_zunmbr
                    call la_zunmbr('Q','L','C',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_zunmbr = real(dum(1),KIND=dp)
                    ! compute space needed for la_zungbr
                    call la_zungbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_zungbr = real(dum(1),KIND=dp)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,2*n + lwork_zgebrd)
                    maxwrk = max(maxwrk,2*n + lwork_zunmbr)
                    maxwrk = max(maxwrk,2*n + lwork_zungbr)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = 2*n + max(nrhs,m)
                 end if
                 if (n > m) then
                    minwrk = 2*m + max(nrhs,n)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_zgelqf
                       call la_zgelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_zgelqf = real(dum(1),KIND=dp)
                       ! compute space needed for la_zgebrd
                       call la_zgebrd(m,m,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_zgebrd = real(dum(1),KIND=dp)
                       ! compute space needed for la_zunmbr
                       call la_zunmbr('Q','L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_zunmbr = real(dum(1),KIND=dp)
                       ! compute space needed for la_zungbr
                       call la_zungbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_zungbr = real(dum(1),KIND=dp)
                       ! compute space needed for la_zunmlq
                       call la_zunmlq('L','C',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_zunmlq = real(dum(1),KIND=dp)
                       ! compute total workspace needed
                       maxwrk = m + lwork_zgelqf
                       maxwrk = max(maxwrk,3*m + m*m + lwork_zgebrd)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_zunmbr)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_zungbr)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_zunmlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_zgebrd
                       call la_zgebrd(m,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_zgebrd = real(dum(1),KIND=dp)
                       ! compute space needed for la_zunmbr
                       call la_zunmbr('Q','L','C',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_zunmbr = real(dum(1),KIND=dp)
                       ! compute space needed for la_zungbr
                       call la_zungbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_zungbr = real(dum(1),KIND=dp)
                       maxwrk = 2*m + lwork_zgebrd
                       maxwrk = max(maxwrk,2*m + lwork_zunmbr)
                       maxwrk = max(maxwrk,2*m + lwork_zungbr)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('ZGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_dlamch('P')
           sfmin = la_dlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_zlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_dlaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (cworkspace: need 2*n, prefer n+n*nb)
                 ! (rworkspace: none)
                 call la_zgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (cworkspace: need n+nrhs, prefer n+nrhs*nb)
                 ! (rworkspace: none)
                 call la_zunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_zlaset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
              end if
              ie = 1
              itauq = 1
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              ! (rworkspace: need n)
              call la_zgebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              ! (rworkspace: none)
              call la_zunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (cworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zungbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_zbdsqr('U',n,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_zdrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_zlaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_zgemm('C','N',n,nrhs,n,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_zlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_zgemm('C','N',n,bl,n,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_zlacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_zgemv('C',n,n,cone,a,lda,b,1,czero,work,1)
                 call la_zcopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 3*m + m*m + max(m,nrhs,n - 2*m)) then
              ! underdetermined case, m much less than n
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= 3*m + m*lda + max(m,nrhs,n - 2*m)) ldwork = lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (cworkspace: need 2*m, prefer m+m*nb)
              ! (rworkspace: none)
              call la_zgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_zlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_zlaset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              ie = 1
              itauq = il + ldwork*m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (cworkspace: need m*m+4*m, prefer m*m+3*m+2*m*nb)
              ! (rworkspace: need m)
              call la_zgebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (cworkspace: need m*m+3*m+nrhs, prefer m*m+3*m+nrhs*nb)
              ! (rworkspace: none)
              call la_zunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (cworkspace: need m*m+4*m-1, prefer m*m+3*m+(m-1)*nb)
              ! (rworkspace: none)
              call la_zungbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration, computing right singular
              ! vectors of l in work(il) and multiplying b by transpose of
              ! left singular vectors
              ! (cworkspace: need m*m)
              ! (rworkspace: need bdspac)
              call la_zbdsqr('U',m,m,0,nrhs,s,rwork(ie),work(il),ldwork,a,lda, &
                        b,ldb,rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_zdrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_zlaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              iwork = il + m*ldwork
              ! multiply b by right singular vectors of l in work(il)
              ! (cworkspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_zgemm('C','N',m,nrhs,m,cone,work(il),ldwork,b,ldb,czero, &
                           work(iwork),ldb)
                 call la_zlacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_zgemm('C','N',m,bl,m,cone,work(il),ldwork,b(1,i), &
                              ldb,czero,work(iwork),m)
                    call la_zlacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_zgemv('C',m,m,cone,work(il),ldwork,b(1,1),1,czero,work( &
                            iwork),1)
                 call la_zcopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_zlaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (cworkspace: need m+nrhs, prefer m+nhrs*nb)
              ! (rworkspace: none)
              call la_zunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = 1
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (cworkspace: need 3*m, prefer 2*m+(m+n)*nb)
              ! (rworkspace: need n)
              call la_zgebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              ! (rworkspace: none)
              call la_zunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (cworkspace: need 3*m, prefer 2*m+m*nb)
              ! (rworkspace: none)
              call la_zungbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_zbdsqr('L',m,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_zdrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_zlaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_zgemm('C','N',n,nrhs,m,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_zlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_zgemm('C','N',n,bl,m,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_zlacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_zgemv('C',m,n,cone,a,lda,b,1,czero,work,1)
                 call la_zcopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_zlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_zlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_dlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_zlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_zlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_zgelss
#ifdef LA_WITH_XDP
     !> YGELSS: computes the minimum norm solution to a complex linear
     !> least squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_ygelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_xdp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(xdp),intent(in) :: rcond
           ! Array Arguments
           real(xdp),intent(out) :: rwork(*),s(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bl,chunk,i,iascl,ibscl,ie,il,irwork,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_ygeqrf,lwork_yunmqr,lwork_ygebrd,lwork_yunmbr,lwork_yungbr, &
                     lwork_yunmlq,lwork_ygelqf
           real(xdp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           complex(xdp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace refers
             ! to real workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'YGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_ygeqrf
                    call la_ygeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_ygeqrf = real(dum(1),KIND=xdp)
                    ! compute space needed for la_yunmqr
                    call la_yunmqr('L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_yunmqr = real(dum(1),KIND=xdp)
                    mm = n
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'YGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'YUNMQR','LC',m,nrhs,n,- &
                              1))
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute space needed for la_ygebrd
                    call la_ygebrd(mm,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                    lwork_ygebrd = real(dum(1),KIND=xdp)
                    ! compute space needed for la_yunmbr
                    call la_yunmbr('Q','L','C',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_yunmbr = real(dum(1),KIND=xdp)
                    ! compute space needed for la_yungbr
                    call la_yungbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_yungbr = real(dum(1),KIND=xdp)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,2*n + lwork_ygebrd)
                    maxwrk = max(maxwrk,2*n + lwork_yunmbr)
                    maxwrk = max(maxwrk,2*n + lwork_yungbr)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = 2*n + max(nrhs,m)
                 end if
                 if (n > m) then
                    minwrk = 2*m + max(nrhs,n)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_ygelqf
                       call la_ygelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_ygelqf = real(dum(1),KIND=xdp)
                       ! compute space needed for la_ygebrd
                       call la_ygebrd(m,m,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_ygebrd = real(dum(1),KIND=xdp)
                       ! compute space needed for la_yunmbr
                       call la_yunmbr('Q','L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_yunmbr = real(dum(1),KIND=xdp)
                       ! compute space needed for la_yungbr
                       call la_yungbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_yungbr = real(dum(1),KIND=xdp)
                       ! compute space needed for la_yunmlq
                       call la_yunmlq('L','C',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_yunmlq = real(dum(1),KIND=xdp)
                       ! compute total workspace needed
                       maxwrk = m + lwork_ygelqf
                       maxwrk = max(maxwrk,3*m + m*m + lwork_ygebrd)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_yunmbr)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_yungbr)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_yunmlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_ygebrd
                       call la_ygebrd(m,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_ygebrd = real(dum(1),KIND=xdp)
                       ! compute space needed for la_yunmbr
                       call la_yunmbr('Q','L','C',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_yunmbr = real(dum(1),KIND=xdp)
                       ! compute space needed for la_yungbr
                       call la_yungbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_yungbr = real(dum(1),KIND=xdp)
                       maxwrk = 2*m + lwork_ygebrd
                       maxwrk = max(maxwrk,2*m + lwork_yunmbr)
                       maxwrk = max(maxwrk,2*m + lwork_yungbr)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('YGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_xlamch('P')
           sfmin = la_xlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_ylange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_ylaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_xlaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_ylange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (cworkspace: need 2*n, prefer n+n*nb)
                 ! (rworkspace: none)
                 call la_ygeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (cworkspace: need n+nrhs, prefer n+nrhs*nb)
                 ! (rworkspace: none)
                 call la_yunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_ylaset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
              end if
              ie = 1
              itauq = 1
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              ! (rworkspace: need n)
              call la_ygebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              ! (rworkspace: none)
              call la_yunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (cworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              ! (rworkspace: none)
              call la_yungbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_ybdsqr('U',n,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_yxrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_ylaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_ygemm('C','N',n,nrhs,n,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_ylacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_ygemm('C','N',n,bl,n,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_ylacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_ygemv('C',n,n,cone,a,lda,b,1,czero,work,1)
                 call la_ycopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 3*m + m*m + max(m,nrhs,n - 2*m)) then
              ! underdetermined case, m much less than n
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= 3*m + m*lda + max(m,nrhs,n - 2*m)) ldwork = lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (cworkspace: need 2*m, prefer m+m*nb)
              ! (rworkspace: none)
              call la_ygelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_ylacpy('L',m,m,a,lda,work(il),ldwork)
              call la_ylaset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              ie = 1
              itauq = il + ldwork*m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (cworkspace: need m*m+4*m, prefer m*m+3*m+2*m*nb)
              ! (rworkspace: need m)
              call la_ygebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (cworkspace: need m*m+3*m+nrhs, prefer m*m+3*m+nrhs*nb)
              ! (rworkspace: none)
              call la_yunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (cworkspace: need m*m+4*m-1, prefer m*m+3*m+(m-1)*nb)
              ! (rworkspace: none)
              call la_yungbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration, computing right singular
              ! vectors of l in work(il) and multiplying b by transpose of
              ! left singular vectors
              ! (cworkspace: need m*m)
              ! (rworkspace: need bdspac)
              call la_ybdsqr('U',m,m,0,nrhs,s,rwork(ie),work(il),ldwork,a,lda, &
                        b,ldb,rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_yxrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_ylaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              iwork = il + m*ldwork
              ! multiply b by right singular vectors of l in work(il)
              ! (cworkspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_ygemm('C','N',m,nrhs,m,cone,work(il),ldwork,b,ldb,czero, &
                           work(iwork),ldb)
                 call la_ylacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_ygemm('C','N',m,bl,m,cone,work(il),ldwork,b(1,i), &
                              ldb,czero,work(iwork),m)
                    call la_ylacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_ygemv('C',m,m,cone,work(il),ldwork,b(1,1),1,czero,work( &
                            iwork),1)
                 call la_ycopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_ylaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (cworkspace: need m+nrhs, prefer m+nhrs*nb)
              ! (rworkspace: none)
              call la_yunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = 1
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (cworkspace: need 3*m, prefer 2*m+(m+n)*nb)
              ! (rworkspace: need n)
              call la_ygebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              ! (rworkspace: none)
              call la_yunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (cworkspace: need 3*m, prefer 2*m+m*nb)
              ! (rworkspace: none)
              call la_yungbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_ybdsqr('L',m,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_yxrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_ylaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_ygemm('C','N',n,nrhs,m,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_ylacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_ygemm('C','N',n,bl,m,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_ylacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_ygemv('C',m,n,cone,a,lda,b,1,czero,work,1)
                 call la_ycopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_ylascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_ylascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_xlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_ylascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_ylascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_ygelss
#endif
#ifdef LA_WITH_QP
     !> WGELSS: computes the minimum norm solution to a complex linear
     !> least squares problem:
     !> Minimize 2-norm(| b - A*x |).
     !> using the singular value decomposition (SVD) of A. A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
     !> X.
     !> The effective rank of A is determined by treating as zero those
     !> singular values which are less than RCOND times the largest singular
     !> value.

     subroutine la_wgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(qp),intent(in) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: rwork(*),s(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: bl,chunk,i,iascl,ibscl,ie,il,irwork,itau,itaup,itauq,iwork, &
                     ldwork,maxmn,maxwrk,minmn,minwrk,mm,mnthr
           integer(ilp) :: lwork_wgeqrf,lwork_wunmqr,lwork_wgebrd,lwork_wunmbr,lwork_wungbr, &
                     lwork_wunmlq,lwork_wgelqf
           real(qp) :: anrm,bignum,bnrm,eps,sfmin,smlnum,thr
           ! Local Arrays
           complex(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           minmn = min(m,n)
           maxmn = max(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,maxmn)) then
              info = -7
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace refers
             ! to real workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              minwrk = 1
              maxwrk = 1
              if (minmn > 0) then
                 mm = m
                 mnthr = la_ilaenv(6,'WGELSS',' ',m,n,nrhs,-1)
                 if (m >= n .and. m >= mnthr) then
                    ! path 1a - overdetermined, with many more rows than
                              ! columns
                    ! compute space needed for la_wgeqrf
                    call la_wgeqrf(m,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_wgeqrf = real(dum(1),KIND=qp)
                    ! compute space needed for la_wunmqr
                    call la_wunmqr('L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum(1),-1, &
                              info)
                    lwork_wunmqr = real(dum(1),KIND=qp)
                    mm = n
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'WGEQRF',' ',m,n,-1,-1))

                    maxwrk = max(maxwrk,n + nrhs*la_ilaenv(1,'WUNMQR','LC',m,nrhs,n,- &
                              1))
                 end if
                 if (m >= n) then
                    ! path 1 - overdetermined or exactly determined
                    ! compute space needed for la_wgebrd
                    call la_wgebrd(mm,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                    lwork_wgebrd = real(dum(1),KIND=qp)
                    ! compute space needed for la_wunmbr
                    call la_wunmbr('Q','L','C',mm,nrhs,n,a,lda,dum(1),b,ldb,dum(1), &
                               -1,info)
                    lwork_wunmbr = real(dum(1),KIND=qp)
                    ! compute space needed for la_wungbr
                    call la_wungbr('P',n,n,n,a,lda,dum(1),dum(1),-1,info)
                    lwork_wungbr = real(dum(1),KIND=qp)
                    ! compute total workspace needed
                    maxwrk = max(maxwrk,2*n + lwork_wgebrd)
                    maxwrk = max(maxwrk,2*n + lwork_wunmbr)
                    maxwrk = max(maxwrk,2*n + lwork_wungbr)
                    maxwrk = max(maxwrk,n*nrhs)
                    minwrk = 2*n + max(nrhs,m)
                 end if
                 if (n > m) then
                    minwrk = 2*m + max(nrhs,n)
                    if (n >= mnthr) then
                       ! path 2a - underdetermined, with many more columns
                       ! than rows
                       ! compute space needed for la_wgelqf
                       call la_wgelqf(m,n,a,lda,dum(1),dum(1),-1,info)
                       lwork_wgelqf = real(dum(1),KIND=qp)
                       ! compute space needed for la_wgebrd
                       call la_wgebrd(m,m,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_wgebrd = real(dum(1),KIND=qp)
                       ! compute space needed for la_wunmbr
                       call la_wunmbr('Q','L','C',m,nrhs,n,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_wunmbr = real(dum(1),KIND=qp)
                       ! compute space needed for la_wungbr
                       call la_wungbr('P',m,m,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_wungbr = real(dum(1),KIND=qp)
                       ! compute space needed for la_wunmlq
                       call la_wunmlq('L','C',n,nrhs,m,a,lda,dum(1),b,ldb,dum(1),- &
                                 1,info)
                       lwork_wunmlq = real(dum(1),KIND=qp)
                       ! compute total workspace needed
                       maxwrk = m + lwork_wgelqf
                       maxwrk = max(maxwrk,3*m + m*m + lwork_wgebrd)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_wunmbr)
                       maxwrk = max(maxwrk,3*m + m*m + lwork_wungbr)
                       if (nrhs > 1) then
                          maxwrk = max(maxwrk,m*m + m + m*nrhs)
                       else
                          maxwrk = max(maxwrk,m*m + 2*m)
                       end if
                       maxwrk = max(maxwrk,m + lwork_wunmlq)
                    else
                       ! path 2 - underdetermined
                       ! compute space needed for la_wgebrd
                       call la_wgebrd(m,n,a,lda,s,s,dum(1),dum(1),dum(1),-1,info)

                       lwork_wgebrd = real(dum(1),KIND=qp)
                       ! compute space needed for la_wunmbr
                       call la_wunmbr('Q','L','C',m,nrhs,m,a,lda,dum(1),b,ldb,dum( &
                                 1),-1,info)
                       lwork_wunmbr = real(dum(1),KIND=qp)
                       ! compute space needed for la_wungbr
                       call la_wungbr('P',m,n,m,a,lda,dum(1),dum(1),-1,info)
                       lwork_wungbr = real(dum(1),KIND=qp)
                       maxwrk = 2*m + lwork_wgebrd
                       maxwrk = max(maxwrk,2*m + lwork_wunmbr)
                       maxwrk = max(maxwrk,2*m + lwork_wungbr)
                       maxwrk = max(maxwrk,n*nrhs)
                    end if
                 end if
                 maxwrk = max(minwrk,maxwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -12
           end if
           if (info /= 0) then
              call la_xerbla('WGELSS',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           eps = la_qlamch('P')
           sfmin = la_qlamch('S')
           smlnum = sfmin/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_wlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              call la_qlaset('F',minmn,1,zero,zero,s,minmn)
              rank = 0
              go to 70
           end if
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! overdetermined case
           if (m >= n) then
              ! path 1 - overdetermined or exactly determined
              mm = m
              if (m >= mnthr) then
                 ! path 1a - overdetermined, with many more rows than columns
                 mm = n
                 itau = 1
                 iwork = itau + n
                 ! compute a=q*r
                 ! (cworkspace: need 2*n, prefer n+n*nb)
                 ! (rworkspace: none)
                 call la_wgeqrf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1, &
                           info)
                 ! multiply b by transpose(q)
                 ! (cworkspace: need n+nrhs, prefer n+nrhs*nb)
                 ! (rworkspace: none)
                 call la_wunmqr('L','C',m,nrhs,n,a,lda,work(itau),b,ldb,work( &
                           iwork),lwork - iwork + 1,info)
                 ! zero out below r
                 if (n > 1) call la_wlaset('L',n - 1,n - 1,czero,czero,a(2,1),lda)
              end if
              ie = 1
              itauq = 1
              itaup = itauq + n
              iwork = itaup + n
              ! bidiagonalize r in a
              ! (cworkspace: need 2*n+mm, prefer 2*n+(mm+n)*nb)
              ! (rworkspace: need n)
              call la_wgebrd(mm,n,a,lda,s,rwork(ie),work(itauq),work(itaup), &
                        work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of r
              ! (cworkspace: need 2*n+nrhs, prefer 2*n+nrhs*nb)
              ! (rworkspace: none)
              call la_wunmbr('Q','L','C',mm,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in a
              ! (cworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wungbr('P',n,n,n,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + n
              ! perform bidiagonal qr iteration
                ! multiply b by transpose of left singular vectors
                ! compute right singular vectors in a
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_wbdsqr('U',n,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,n
                 if (s(i) > thr) then
                    call la_wqrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_wlaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_wgemm('C','N',n,nrhs,n,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_wlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_wgemm('C','N',n,bl,n,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_wlacpy('G',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_wgemv('C',n,n,cone,a,lda,b,1,czero,work,1)
                 call la_wcopy(n,work,1,b,1)
              end if
           else if (n >= mnthr .and. lwork >= 3*m + m*m + max(m,nrhs,n - 2*m)) then
              ! underdetermined case, m much less than n
              ! path 2a - underdetermined, with many more columns than rows
              ! and sufficient workspace for an efficient algorithm
              ldwork = m
              if (lwork >= 3*m + m*lda + max(m,nrhs,n - 2*m)) ldwork = lda
              itau = 1
              iwork = m + 1
              ! compute a=l*q
              ! (cworkspace: need 2*m, prefer m+m*nb)
              ! (rworkspace: none)
              call la_wgelqf(m,n,a,lda,work(itau),work(iwork),lwork - iwork + 1,info)

              il = iwork
              ! copy l to work(il), zeroing out above it
              call la_wlacpy('L',m,m,a,lda,work(il),ldwork)
              call la_wlaset('U',m - 1,m - 1,czero,czero,work(il + ldwork),ldwork)
              ie = 1
              itauq = il + ldwork*m
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize l in work(il)
              ! (cworkspace: need m*m+4*m, prefer m*m+3*m+2*m*nb)
              ! (rworkspace: need m)
              call la_wgebrd(m,m,work(il),ldwork,s,rwork(ie),work(itauq),work( &
                        itaup),work(iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors of l
              ! (cworkspace: need m*m+3*m+nrhs, prefer m*m+3*m+nrhs*nb)
              ! (rworkspace: none)
              call la_wunmbr('Q','L','C',m,nrhs,m,work(il),ldwork,work(itauq),b, &
                        ldb,work(iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors of r in work(il)
              ! (cworkspace: need m*m+4*m-1, prefer m*m+3*m+(m-1)*nb)
              ! (rworkspace: none)
              call la_wungbr('P',m,m,m,work(il),ldwork,work(itaup),work(iwork), &
                        lwork - iwork + 1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration, computing right singular
              ! vectors of l in work(il) and multiplying b by transpose of
              ! left singular vectors
              ! (cworkspace: need m*m)
              ! (rworkspace: need bdspac)
              call la_wbdsqr('U',m,m,0,nrhs,s,rwork(ie),work(il),ldwork,a,lda, &
                        b,ldb,rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_wqrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_wlaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              iwork = il + m*ldwork
              ! multiply b by right singular vectors of l in work(il)
              ! (cworkspace: need m*m+2*m, prefer m*m+m+m*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs + iwork - 1 .and. nrhs > 1) then
                 call la_wgemm('C','N',m,nrhs,m,cone,work(il),ldwork,b,ldb,czero, &
                           work(iwork),ldb)
                 call la_wlacpy('G',m,nrhs,work(iwork),ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = (lwork - iwork + 1)/m
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_wgemm('C','N',m,bl,m,cone,work(il),ldwork,b(1,i), &
                              ldb,czero,work(iwork),m)
                    call la_wlacpy('G',m,bl,work(iwork),m,b(1,i),ldb)
                 end do
              else
                 call la_wgemv('C',m,m,cone,work(il),ldwork,b(1,1),1,czero,work( &
                            iwork),1)
                 call la_wcopy(m,work(iwork),1,b(1,1),1)
              end if
              ! zero out below first m rows of b
              call la_wlaset('F',n - m,nrhs,czero,czero,b(m + 1,1),ldb)
              iwork = itau + m
              ! multiply transpose(q) by b
              ! (cworkspace: need m+nrhs, prefer m+nhrs*nb)
              ! (rworkspace: none)
              call la_wunmlq('L','C',n,nrhs,m,a,lda,work(itau),b,ldb,work(iwork) &
                        ,lwork - iwork + 1,info)
           else
              ! path 2 - remaining underdetermined cases
              ie = 1
              itauq = 1
              itaup = itauq + m
              iwork = itaup + m
              ! bidiagonalize a
              ! (cworkspace: need 3*m, prefer 2*m+(m+n)*nb)
              ! (rworkspace: need n)
              call la_wgebrd(m,n,a,lda,s,rwork(ie),work(itauq),work(itaup),work( &
                         iwork),lwork - iwork + 1,info)
              ! multiply b by transpose of left bidiagonalizing vectors
              ! (cworkspace: need 2*m+nrhs, prefer 2*m+nrhs*nb)
              ! (rworkspace: none)
              call la_wunmbr('Q','L','C',m,nrhs,n,a,lda,work(itauq),b,ldb,work( &
                        iwork),lwork - iwork + 1,info)
              ! generate right bidiagonalizing vectors in a
              ! (cworkspace: need 3*m, prefer 2*m+m*nb)
              ! (rworkspace: none)
              call la_wungbr('P',m,n,m,a,lda,work(itaup),work(iwork),lwork - iwork + &
                        1,info)
              irwork = ie + m
              ! perform bidiagonal qr iteration,
                 ! computing right singular vectors of a in a and
                 ! multiplying b by transpose of left singular vectors
              ! (cworkspace: none)
              ! (rworkspace: need bdspac)
              call la_wbdsqr('L',m,n,0,nrhs,s,rwork(ie),a,lda,dum,1,b,ldb, &
                        rwork(irwork),info)
              if (info /= 0) go to 70
              ! multiply b by reciprocals of singular values
              thr = max(rcond*s(1),sfmin)
              if (rcond < zero) thr = max(eps*s(1),sfmin)
              rank = 0
              do i = 1,m
                 if (s(i) > thr) then
                    call la_wqrscl(nrhs,s(i),b(i,1),ldb)
                    rank = rank + 1
                 else
                    call la_wlaset('F',1,nrhs,czero,czero,b(i,1),ldb)
                 end if
              end do
              ! multiply b by right singular vectors of a
              ! (cworkspace: need n, prefer n*nrhs)
              ! (rworkspace: none)
              if (lwork >= ldb*nrhs .and. nrhs > 1) then
                 call la_wgemm('C','N',n,nrhs,m,cone,a,lda,b,ldb,czero,work,ldb)

                 call la_wlacpy('G',n,nrhs,work,ldb,b,ldb)
              else if (nrhs > 1) then
                 chunk = lwork/n
                 do i = 1,nrhs,chunk
                    bl = min(nrhs - i + 1,chunk)
                    call la_wgemm('C','N',n,bl,m,cone,a,lda,b(1,i),ldb,czero, &
                              work,n)
                    call la_wlacpy('F',n,bl,work,n,b(1,i),ldb)
                 end do
              else
                 call la_wgemv('C',m,n,cone,a,lda,b,1,czero,work,1)
                 call la_wcopy(n,work,1,b,1)
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
              call la_wlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,smlnum,anrm,minmn,1,s,minmn,info)
           else if (iascl == 2) then
              call la_wlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_qlascl('G',0,0,bignum,anrm,minmn,1,s,minmn,info)
           end if
           if (ibscl == 1) then
              call la_wlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_wlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = maxwrk
           return
     end subroutine la_wgelss
#endif

     !> CGELSY: computes the minimum-norm solution to a complex linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by unitary transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**H [ inv(T11)*Q1**H*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.

     subroutine la_cgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(sp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkopt,mn,nb,nb1,nb2,nb3, &
                     nb4
           real(sp) :: anrm,bignum,bnrm,smax,smaxpr,smin,sminpr,smlnum,wsize
           complex(sp) :: c1,c2,s1,s2
           ! Intrinsic Functions
           intrinsic :: abs,max,min,real,cmplx
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           nb1 = la_ilaenv(1,'CGEQRF',' ',m,n,-1,-1)
           nb2 = la_ilaenv(1,'CGERQF',' ',m,n,-1,-1)
           nb3 = la_ilaenv(1,'CUNMQR',' ',m,n,nrhs,-1)
           nb4 = la_ilaenv(1,'CUNMRQ',' ',m,n,nrhs,-1)
           nb = max(nb1,nb2,nb3,nb4)
           lwkopt = max(1,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
           work(1) = cmplx(lwkopt,KIND=sp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           else if (lwork < (mn + max(2*mn,n + 1,mn + nrhs)) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('CGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_slamch('S')/la_slamch('P')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_clange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_claset('F',max(m,n),nrhs,czero,czero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_clange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_cgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,rwork,info)

           wsize = mn + real(work(mn + 1),KIND=sp)
           ! complex workspace: mn+nb*(n+1). real workspace 2*n.
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = cone
           work(ismax) = cone
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_claset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_claic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_claic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! complex workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_ctzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! complex workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
           call la_cunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,mn,a,lda,work(1),b, &
                     ldb,work(2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + real(work(2*mn + 1),KIND=sp))
           ! complex workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_ctrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,cone,a, &
                     lda,b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = czero
              end do
           end do
           ! b(1:n,1:nrhs) := y**h * b(1:n,1:nrhs)
           if (rank < n) then
              call la_cunmrz('LEFT','CONJUGATE TRANSPOSE',n,nrhs,rank,n - rank,a,lda, &
                        work(mn + 1),b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! complex workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_ccopy(n,work(1),1,b(1,j),1)
           end do
           ! complex workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_clascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_clascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_clascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_clascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_clascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_clascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = cmplx(lwkopt,KIND=sp)
           return
     end subroutine la_cgelsy
     !> ZGELSY: computes the minimum-norm solution to a complex linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by unitary transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**H [ inv(T11)*Q1**H*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.

     subroutine la_zgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(dp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkopt,mn,nb,nb1,nb2,nb3, &
                     nb4
           real(dp) :: anrm,bignum,bnrm,smax,smaxpr,smin,sminpr,smlnum,wsize
           complex(dp) :: c1,c2,s1,s2
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           nb1 = la_ilaenv(1,'ZGEQRF',' ',m,n,-1,-1)
           nb2 = la_ilaenv(1,'ZGERQF',' ',m,n,-1,-1)
           nb3 = la_ilaenv(1,'ZUNMQR',' ',m,n,nrhs,-1)
           nb4 = la_ilaenv(1,'ZUNMRQ',' ',m,n,nrhs,-1)
           nb = max(nb1,nb2,nb3,nb4)
           lwkopt = max(1,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
           work(1) = cmplx(lwkopt,KIND=dp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           else if (lwork < (mn + max(2*mn,n + 1,mn + nrhs)) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('ZGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_dlamch('S')/la_dlamch('P')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_zlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_zlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_zlange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_zgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,rwork,info)

           wsize = mn + real(work(mn + 1),KIND=dp)
           ! complex workspace: mn+nb*(n+1). real workspace 2*n.
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = cone
           work(ismax) = cone
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_zlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_zlaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_zlaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! complex workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_ztzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! complex workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
           call la_zunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,mn,a,lda,work(1),b, &
                     ldb,work(2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + real(work(2*mn + 1),KIND=dp))
           ! complex workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_ztrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,cone,a, &
                     lda,b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = czero
              end do
           end do
           ! b(1:n,1:nrhs) := y**h * b(1:n,1:nrhs)
           if (rank < n) then
              call la_zunmrz('LEFT','CONJUGATE TRANSPOSE',n,nrhs,rank,n - rank,a,lda, &
                        work(mn + 1),b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! complex workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_zcopy(n,work(1),1,b(1,j),1)
           end do
           ! complex workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_zlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_zlascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_zlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_zlascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_zlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_zlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = cmplx(lwkopt,KIND=dp)
           return
     end subroutine la_zgelsy
#ifdef LA_WITH_XDP
     !> YGELSY: computes the minimum-norm solution to a complex linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by unitary transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**H [ inv(T11)*Q1**H*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.

     subroutine la_ygelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_xdp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(xdp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkopt,mn,nb,nb1,nb2,nb3, &
                     nb4
           real(xdp) :: anrm,bignum,bnrm,smax,smaxpr,smin,sminpr,smlnum,wsize
           complex(xdp) :: c1,c2,s1,s2
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           nb1 = la_ilaenv(1,'YGEQRF',' ',m,n,-1,-1)
           nb2 = la_ilaenv(1,'YGERQF',' ',m,n,-1,-1)
           nb3 = la_ilaenv(1,'YUNMQR',' ',m,n,nrhs,-1)
           nb4 = la_ilaenv(1,'YUNMRQ',' ',m,n,nrhs,-1)
           nb = max(nb1,nb2,nb3,nb4)
           lwkopt = max(1,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
           work(1) = cmplx(lwkopt,KIND=xdp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           else if (lwork < (mn + max(2*mn,n + 1,mn + nrhs)) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('YGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_xlamch('S')/la_xlamch('P')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_ylange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_ylaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_ylange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_ygeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,rwork,info)

           wsize = mn + real(work(mn + 1),KIND=xdp)
           ! complex workspace: mn+nb*(n+1). real workspace 2*n.
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = cone
           work(ismax) = cone
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_ylaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_ylaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_ylaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! complex workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_ytzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! complex workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
           call la_yunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,mn,a,lda,work(1),b, &
                     ldb,work(2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + real(work(2*mn + 1),KIND=xdp))
           ! complex workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_ytrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,cone,a, &
                     lda,b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = czero
              end do
           end do
           ! b(1:n,1:nrhs) := y**h * b(1:n,1:nrhs)
           if (rank < n) then
              call la_yunmrz('LEFT','CONJUGATE TRANSPOSE',n,nrhs,rank,n - rank,a,lda, &
                        work(mn + 1),b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! complex workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_ycopy(n,work(1),1,b(1,j),1)
           end do
           ! complex workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_ylascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_ylascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_ylascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_ylascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_ylascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_ylascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = cmplx(lwkopt,KIND=xdp)
           return
     end subroutine la_ygelsy
#endif
#ifdef LA_WITH_QP
     !> WGELSY: computes the minimum-norm solution to a complex linear least
     !> squares problem:
     !> minimize || A * X - B ||
     !> using a complete orthogonal factorization of A.  A is an M-by-N
     !> matrix which may be rank-deficient.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.
     !> The routine first computes a QR factorization with column pivoting:
     !> A * P = Q * [ R11 R12 ]
     !> [  0  R22 ]
     !> with R11 defined as the largest leading submatrix whose estimated
     !> condition number is less than 1/RCOND.  The order of R11, RANK,
     !> is the effective rank of A.
     !> Then, R22 is considered to be negligible, and R12 is annihilated
     !> by unitary transformations from the right, arriving at the
     !> complete orthogonal factorization:
     !> A * P = Q * [ T11 0 ] * Z
     !> [  0  0 ]
     !> The minimum-norm solution is then
     !> X = P * Z**H [ inv(T11)*Q1**H*B ]
     !> [        0         ]
     !> where Q1 consists of the first RANK columns of Q.
     !> This routine is basically identical to the original xGELSX except
     !> three differences:
     !> o The permutation of matrix B (the right hand side) is faster and
     !> more simple.
     !> o The call to the subroutine xGEQPF has been substituted by the
     !> the call to the subroutine xGEQP3. This subroutine is a Blas-3
     !> version of the QR factorization with column pivoting.
     !> o Matrix B (the right hand side) is updated with Blas-3.

     subroutine la_wgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,rwork, &
               info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,rank
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           real(qp),intent(in) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: jpvt(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: imax = 1
           integer(ilp),parameter :: imin = 2

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iascl,ibscl,ismax,ismin,j,lwkopt,mn,nb,nb1,nb2,nb3, &
                     nb4
           real(qp) :: anrm,bignum,bnrm,smax,smaxpr,smin,sminpr,smlnum,wsize
           complex(qp) :: c1,c2,s1,s2
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,max,min
           ! Executable Statements
           mn = min(m,n)
           ismin = mn + 1
           ismax = 2*mn + 1
           ! test the input arguments.
           info = 0
           nb1 = la_ilaenv(1,'WGEQRF',' ',m,n,-1,-1)
           nb2 = la_ilaenv(1,'WGERQF',' ',m,n,-1,-1)
           nb3 = la_ilaenv(1,'WUNMQR',' ',m,n,nrhs,-1)
           nb4 = la_ilaenv(1,'WUNMRQ',' ',m,n,nrhs,-1)
           nb = max(nb1,nb2,nb3,nb4)
           lwkopt = max(1,mn + 2*n + nb*(n + 1),2*mn + nb*nrhs)
           work(1) = cmplx(lwkopt,KIND=qp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m,n)) then
              info = -7
           else if (lwork < (mn + max(2*mn,n + 1,mn + nrhs)) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('WGELSY',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
              rank = 0
              return
           end if
           ! get machine parameters
           smlnum = la_qlamch('S')/la_qlamch('P')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! scale a, b if max entries outside range [smlnum,bignum]
           anrm = la_wlange('M',m,n,a,lda,rwork)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_wlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              rank = 0
              go to 70
           end if
           bnrm = la_wlange('M',m,nrhs,b,ldb,rwork)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,bnrm,smlnum,m,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,bnrm,bignum,m,nrhs,b,ldb,info)
              ibscl = 2
           end if
           ! compute qr factorization with column pivoting of a:
              ! a * p = q * r
           call la_wgeqp3(m,n,a,lda,jpvt,work(1),work(mn + 1),lwork - mn,rwork,info)

           wsize = mn + real(work(mn + 1),KIND=qp)
           ! complex workspace: mn+nb*(n+1). real workspace 2*n.
           ! details of householder rotations stored in work(1:mn).
           ! determine rank using incremental condition estimation
           work(ismin) = cone
           work(ismax) = cone
           smax = abs(a(1,1))
           smin = smax
           if (abs(a(1,1)) == zero) then
              rank = 0
              call la_wlaset('F',max(m,n),nrhs,czero,czero,b,ldb)
              go to 70
           else
              rank = 1
           end if
           10 continue
           if (rank < mn) then
              i = rank + 1
              call la_wlaic1(imin,rank,work(ismin),smin,a(1,i),a(i,i),sminpr, &
                        s1,c1)
              call la_wlaic1(imax,rank,work(ismax),smax,a(1,i),a(i,i),smaxpr, &
                        s2,c2)
              if (smaxpr*rcond <= sminpr) then
                 do i = 1,rank
                    work(ismin + i - 1) = s1*work(ismin + i - 1)
                    work(ismax + i - 1) = s2*work(ismax + i - 1)
                 end do
                 work(ismin + rank) = c1
                 work(ismax + rank) = c2
                 smin = sminpr
                 smax = smaxpr
                 rank = rank + 1
                 go to 10
              end if
           end if
           ! complex workspace: 3*mn.
           ! logically partition r = [ r11 r12 ]
                                   ! [  0  r22 ]
           ! where r11 = r(1:rank,1:rank)
           ! [r11,r12] = [ t11, 0 ] * y
           if (rank < n) call la_wtzrzf(rank,n,a,lda,work(mn + 1),work(2*mn + 1),lwork - &
                     2*mn,info)
           ! complex workspace: 2*mn.
           ! details of householder rotations stored in work(mn+1:2*mn)
           ! b(1:m,1:nrhs) := q**h * b(1:m,1:nrhs)
           call la_wunmqr('LEFT','CONJUGATE TRANSPOSE',m,nrhs,mn,a,lda,work(1),b, &
                     ldb,work(2*mn + 1),lwork - 2*mn,info)
           wsize = max(wsize,2*mn + real(work(2*mn + 1),KIND=qp))
           ! complex workspace: 2*mn+nb*nrhs.
           ! b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
           call la_wtrsm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',rank,nrhs,cone,a, &
                     lda,b,ldb)
           do j = 1,nrhs
              do i = rank + 1,n
                 b(i,j) = czero
              end do
           end do
           ! b(1:n,1:nrhs) := y**h * b(1:n,1:nrhs)
           if (rank < n) then
              call la_wunmrz('LEFT','CONJUGATE TRANSPOSE',n,nrhs,rank,n - rank,a,lda, &
                        work(mn + 1),b,ldb,work(2*mn + 1),lwork - 2*mn,info)
           end if
           ! complex workspace: 2*mn+nrhs.
           ! b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
           do j = 1,nrhs
              do i = 1,n
                 work(jpvt(i)) = b(i,j)
              end do
              call la_wcopy(n,work(1),1,b(1,j),1)
           end do
           ! complex workspace: n.
           ! undo scaling
           if (iascl == 1) then
              call la_wlascl('G',0,0,anrm,smlnum,n,nrhs,b,ldb,info)
              call la_wlascl('U',0,0,smlnum,anrm,rank,rank,a,lda,info)
           else if (iascl == 2) then
              call la_wlascl('G',0,0,anrm,bignum,n,nrhs,b,ldb,info)
              call la_wlascl('U',0,0,bignum,anrm,rank,rank,a,lda,info)
           end if
           if (ibscl == 1) then
              call la_wlascl('G',0,0,smlnum,bnrm,n,nrhs,b,ldb,info)
           else if (ibscl == 2) then
              call la_wlascl('G',0,0,bignum,bnrm,n,nrhs,b,ldb,info)
           end if
           70 continue
           work(1) = cmplx(lwkopt,KIND=qp)
           return
     end subroutine la_wgelsy
#endif

     !> CGETSLS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_cgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_sp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(sp) :: anrm,bignum,bnrm,smlnum,dum(1)
           complex(sp) :: tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'C')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_cgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_cgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_cgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_cgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_cgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_cgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_cgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_cgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=sp)
           end if
           if (info /= 0) then
             call la_xerbla('CGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=sp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_claset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_slamch('S')/la_slamch('P')
            bignum = one/smlnum
            call la_slabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_clange('M',m,n,a,lda,dum)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_claset('F',maxmn,nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_clange('M',brow,nrhs,b,ldb,dum)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_clascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_clascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_cgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_cgemqr('L','C',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_ctrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_ctrtrs('U','C','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = czero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_cgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_cgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_ctrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_cgemlq('L','C',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_cgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_ctrtrs('L','C','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_clascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_clascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_clascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_clascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=sp)
           return
     end subroutine la_cgetsls
     !> ZGETSLS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_zgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_dp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(dp) :: anrm,bignum,bnrm,smlnum,dum(1)
           complex(dp) :: tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'C')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_zgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_zgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_zgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_zgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_zgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_zgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_zgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_zgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=dp)
           end if
           if (info /= 0) then
             call la_xerbla('ZGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=dp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_zlaset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_dlamch('S')/la_dlamch('P')
            bignum = one/smlnum
            call la_dlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',m,n,a,lda,dum)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_zlaset('F',maxmn,nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_zlange('M',brow,nrhs,b,ldb,dum)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_zlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_zlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_zgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_zgemqr('L','C',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_ztrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_ztrtrs('U','C','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = czero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_zgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_zgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_ztrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_zgemlq('L','C',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_zgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_ztrtrs('L','C','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_zlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_zlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_zlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_zlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=dp)
           return
     end subroutine la_zgetsls
#ifdef LA_WITH_XDP
     !> YGETSLS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_ygetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_xdp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(xdp) :: anrm,bignum,bnrm,smlnum,dum(1)
           complex(xdp) :: tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'C')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_ygeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_ygemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_ygeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_ygemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_ygelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_ygemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_ygelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_ygemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=xdp)
           end if
           if (info /= 0) then
             call la_xerbla('YGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=xdp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_ylaset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_xlamch('S')/la_xlamch('P')
            bignum = one/smlnum
            call la_xlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_ylange('M',m,n,a,lda,dum)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_ylaset('F',maxmn,nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_ylange('M',brow,nrhs,b,ldb,dum)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_ylascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_ylascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_ygeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_ygemqr('L','C',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_ytrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_ytrtrs('U','C','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = czero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_ygemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_ygelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_ytrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_ygemlq('L','C',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_ygemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_ytrtrs('L','C','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_ylascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_ylascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_ylascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_ylascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=xdp)
           return
     end subroutine la_ygetsls
#endif
#ifdef LA_WITH_QP
     !> WGETSLS: solves overdetermined or underdetermined complex linear systems
     !> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     !> factorization of A.  It is assumed that A has full rank.
     !> The following options are provided:
     !> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A*X ||.
     !> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     !> an underdetermined system A * X = B.
     !> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
     !> an undetermined system A**T * X = B.
     !> 4. If TRANS = 'C' and m < n:  find the least squares solution of
     !> an overdetermined system, i.e., solve the least squares problem
     !> minimize || B - A**T * X ||.
     !> Several right hand side vectors b and solution vectors x can be
     !> handled in a single call; they are stored as the columns of the
     !> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     !> matrix X.

     subroutine la_wgetsls(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
        use la_constants_qp,only:zero,one,czero
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,nrhs
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,tran
           integer(ilp) :: i,iascl,ibscl,j,maxmn,brow,scllen,tszo,tszm,lwo,lwm,lw1, &
                     lw2,wsizeo,wsizem,info2
           real(qp) :: anrm,bignum,bnrm,smlnum,dum(1)
           complex(qp) :: tq(5),workq(1)
           ! Intrinsic Functions
           intrinsic :: real,max,min,int
           ! Executable Statements
           ! test the input arguments.
           info = 0
           maxmn = max(m,n)
           tran = la_lsame(trans,'C')
           lquery = (lwork == -1 .or. lwork == -2)
           if (.not. (la_lsame(trans,'N') .or. la_lsame(trans,'C'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m,n)) then
              info = -8
           end if
           if (info == 0) then
           ! determine the optimum and minimum lwork
            if (m >= n) then
              call la_wgeqr(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_wgemqr('L',trans,m,nrhs,n,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_wgeqr(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_wgemqr('L',trans,m,nrhs,n,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            else
              call la_wgelq(m,n,a,lda,tq,-1,workq,-1,info2)
              tszo = int(tq(1),KIND=ilp)
              lwo = int(workq(1),KIND=ilp)
              call la_wgemlq('L',trans,n,nrhs,m,a,lda,tq,tszo,b,ldb,workq,-1, &
                        info2)
              lwo = max(lwo,int(workq(1),KIND=ilp))
              call la_wgelq(m,n,a,lda,tq,-2,workq,-2,info2)
              tszm = int(tq(1),KIND=ilp)
              lwm = int(workq(1),KIND=ilp)
              call la_wgemlq('L',trans,n,nrhs,m,a,lda,tq,tszm,b,ldb,workq,-1, &
                        info2)
              lwm = max(lwm,int(workq(1),KIND=ilp))
              wsizeo = tszo + lwo
              wsizem = tszm + lwm
            end if
            if ((lwork < wsizem) .and. (.not. lquery)) then
               info = -10
            end if
            work(1) = real(wsizeo,KIND=qp)
           end if
           if (info /= 0) then
             call la_xerbla('WGETSLS',-info)
             return
           end if
           if (lquery) then
             if (lwork == -2) work(1) = real(wsizem,KIND=qp)
             return
           end if
           if (lwork < wsizeo) then
             lw1 = tszm
             lw2 = lwm
           else
             lw1 = tszo
             lw2 = lwo
           end if
           ! quick return if possible
           if (min(m,n,nrhs) == 0) then
                call la_wlaset('FULL',max(m,n),nrhs,czero,czero,b,ldb)
                return
           end if
           ! get machine parameters
            smlnum = la_qlamch('S')/la_qlamch('P')
            bignum = one/smlnum
            call la_qlabad(smlnum,bignum)
           ! scale a, b if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',m,n,a,lda,dum)
           iascl = 0
           if (anrm > zero .and. anrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,anrm,smlnum,m,n,a,lda,info)
              iascl = 1
           else if (anrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,anrm,bignum,m,n,a,lda,info)
              iascl = 2
           else if (anrm == zero) then
              ! matrix all zero. return zero solution.
              call la_wlaset('F',maxmn,nrhs,czero,czero,b,ldb)
              go to 50
           end if
           brow = m
           if (tran) then
             brow = n
           end if
           bnrm = la_wlange('M',brow,nrhs,b,ldb,dum)
           ibscl = 0
           if (bnrm > zero .and. bnrm < smlnum) then
              ! scale matrix norm up to smlnum
              call la_wlascl('G',0,0,bnrm,smlnum,brow,nrhs,b,ldb,info)
              ibscl = 1
           else if (bnrm > bignum) then
              ! scale matrix norm down to bignum
              call la_wlascl('G',0,0,bnrm,bignum,brow,nrhs,b,ldb,info)
              ibscl = 2
           end if
           if (m >= n) then
              ! compute qr factorization of a
             call la_wgeqr(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
             if (.not. tran) then
                 ! least-squares problem min || a * x - b ||
                 ! b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
               call la_wgemqr('L','C',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb,work( &
                          1),lw2,info)
                 ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
               call la_wtrtrs('U','N','N',n,nrhs,a,lda,b,ldb,info)
               if (info > 0) then
                 return
               end if
               scllen = n
             else
                 ! overdetermined system of equations a**t * x = b
                 ! b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
                 call la_wtrtrs('U','C','N',n,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(n+1:m,1:nrhs) = czero
                 do j = 1,nrhs
                    do i = n + 1,m
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
                 call la_wgemqr('L','N',m,nrhs,n,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 scllen = m
              end if
           else
              ! compute lq factorization of a
              call la_wgelq(m,n,a,lda,work(lw2 + 1),lw1,work(1),lw2,info)
              ! workspace at least m, optimally m*nb.
              if (.not. tran) then
                 ! underdetermined system of equations a * x = b
                 ! b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
                 call la_wtrtrs('L','N','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 ! b(m+1:n,1:nrhs) = 0
                 do j = 1,nrhs
                    do i = m + 1,n
                       b(i,j) = czero
                    end do
                 end do
                 ! b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
                 call la_wgemlq('L','C',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 scllen = n
              else
                 ! overdetermined system min || a**t * x - b ||
                 ! b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
                 call la_wgemlq('L','N',n,nrhs,m,a,lda,work(lw2 + 1),lw1,b,ldb, &
                           work(1),lw2,info)
                 ! workspace at least nrhs, optimally nrhs*nb
                 ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
                 call la_wtrtrs('L','C','N',m,nrhs,a,lda,b,ldb,info)
                 if (info > 0) then
                    return
                 end if
                 scllen = m
              end if
           end if
           ! undo scaling
           if (iascl == 1) then
             call la_wlascl('G',0,0,anrm,smlnum,scllen,nrhs,b,ldb,info)
           else if (iascl == 2) then
             call la_wlascl('G',0,0,anrm,bignum,scllen,nrhs,b,ldb,info)
           end if
           if (ibscl == 1) then
             call la_wlascl('G',0,0,smlnum,bnrm,scllen,nrhs,b,ldb,info)
           else if (ibscl == 2) then
             call la_wlascl('G',0,0,bignum,bnrm,scllen,nrhs,b,ldb,info)
           end if
           50 continue
           work(1) = real(tszo + lwo,KIND=qp)
           return
     end subroutine la_wgetsls
#endif

end module la_lapack_lsq
