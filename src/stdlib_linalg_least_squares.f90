

module stdlib_linalg_least_squares
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Compute a least squares solution to system Ax=b, i.e. such that the 2-norm abs(b-Ax) is minimized.
     public :: lstsq

     ! NumPy: lstsq(a, b, rcond='warn')
     ! Scipy: lstsq(a, b, cond=None, overwrite_a=False, overwrite_b=False, check_finite=True, lapack_driver=None)
     ! IMSL: Result = IMSL_QRSOL(B, [A] [, AUXQR] [, BASIS] [, /DOUBLE] [, QR] [, PIVOT] [, RESIDUAL] [, TOLERANCE])

     interface lstsq
        module procedure stdlib_linalg_slstsq_one
        module procedure stdlib_linalg_dlstsq_one
        module procedure stdlib_linalg_qlstsq_one
        module procedure stdlib_linalg_clstsq_one
        module procedure stdlib_linalg_zlstsq_one
        module procedure stdlib_linalg_wlstsq_one
        module procedure stdlib_linalg_slstsq_multiple
        module procedure stdlib_linalg_dlstsq_multiple
        module procedure stdlib_linalg_qlstsq_multiple
        module procedure stdlib_linalg_clstsq_multiple
        module procedure stdlib_linalg_zlstsq_multiple
        module procedure stdlib_linalg_wlstsq_multiple
     end interface lstsq


     contains

     ! Workspace needed by gesv
     subroutine sgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp), intent(in) :: m,n,nrhs
         integer(ilp), intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = stdlib_ilaenv(9,'sgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0, ilog2(mnmin/(smlsiz+1))+1)

         ! Real space
         lrwork = 12*mnmin+2*mnmin*smlsiz+8*mnmin*nlvl+mnmin*nrhs+(smlsiz+1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1, 3*mnmin*nlvl+11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine sgesv_space

     ! Workspace needed by gesv
     subroutine dgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp), intent(in) :: m,n,nrhs
         integer(ilp), intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = stdlib_ilaenv(9,'dgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0, ilog2(mnmin/(smlsiz+1))+1)

         ! Real space
         lrwork = 12*mnmin+2*mnmin*smlsiz+8*mnmin*nlvl+mnmin*nrhs+(smlsiz+1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1, 3*mnmin*nlvl+11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine dgesv_space

     ! Workspace needed by gesv
     subroutine qgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp), intent(in) :: m,n,nrhs
         integer(ilp), intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = stdlib_ilaenv(9,'qgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0, ilog2(mnmin/(smlsiz+1))+1)

         ! Real space
         lrwork = 12*mnmin+2*mnmin*smlsiz+8*mnmin*nlvl+mnmin*nrhs+(smlsiz+1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1, 3*mnmin*nlvl+11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine qgesv_space

     ! Workspace needed by gesv
     subroutine cgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp), intent(in) :: m,n,nrhs
         integer(ilp), intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = stdlib_ilaenv(9,'cgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0, ilog2(mnmin/(smlsiz+1))+1)

         ! Real space
         lrwork = 10*mnmin+2*mnmin*smlsiz+8*mnmin*nlvl+3*smlsiz*nrhs+max((smlsiz+1)**2,n*(1+nrhs)+2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1, 3*mnmin*nlvl+11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine cgesv_space

     ! Workspace needed by gesv
     subroutine zgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp), intent(in) :: m,n,nrhs
         integer(ilp), intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = stdlib_ilaenv(9,'zgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0, ilog2(mnmin/(smlsiz+1))+1)

         ! Real space
         lrwork = 10*mnmin+2*mnmin*smlsiz+8*mnmin*nlvl+3*smlsiz*nrhs+max((smlsiz+1)**2,n*(1+nrhs)+2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1, 3*mnmin*nlvl+11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine zgesv_space

     ! Workspace needed by gesv
     subroutine wgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp), intent(in) :: m,n,nrhs
         integer(ilp), intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = stdlib_ilaenv(9,'wgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0, ilog2(mnmin/(smlsiz+1))+1)

         ! Real space
         lrwork = 10*mnmin+2*mnmin*smlsiz+8*mnmin*nlvl+3*smlsiz*nrhs+max((smlsiz+1)**2,n*(1+nrhs)+2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1, 3*mnmin*nlvl+11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine wgesv_space



     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_slstsq_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(sp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),                     intent(in)            :: b(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp), allocatable, target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp), allocatable :: singular(:),rwork(:)
         real(sp), pointer :: xmat(:,:),amat(:,:)
         real(sp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call sgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_slstsq_one


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_dlstsq_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(dp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),                     intent(in)            :: b(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp), allocatable, target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp), allocatable :: singular(:),rwork(:)
         real(dp), pointer :: xmat(:,:),amat(:,:)
         real(dp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call dgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_dlstsq_one


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_qlstsq_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(qp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),                     intent(in)            :: b(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp), allocatable, target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp), allocatable :: singular(:),rwork(:)
         real(qp), pointer :: xmat(:,:),amat(:,:)
         real(qp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call qgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_qlstsq_one


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_clstsq_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),                     intent(in)            :: b(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp), allocatable, target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp), allocatable :: singular(:),rwork(:)
         complex(sp), pointer :: xmat(:,:),amat(:,:)
         complex(sp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call cgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_clstsq_one


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_zlstsq_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),                     intent(in)            :: b(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp), allocatable, target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp), allocatable :: singular(:),rwork(:)
         complex(dp), pointer :: xmat(:,:),amat(:,:)
         complex(dp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call zgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_zlstsq_one


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_wlstsq_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),                     intent(in)            :: b(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp), allocatable, target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp), allocatable :: singular(:),rwork(:)
         complex(qp), pointer :: xmat(:,:),amat(:,:)
         complex(qp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call wgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_wlstsq_one


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_slstsq_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(sp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),                     intent(in)            :: b(:,:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp), allocatable, target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp), allocatable :: singular(:),rwork(:)
         real(sp), pointer :: xmat(:,:),amat(:,:)
         real(sp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call sgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_slstsq_multiple


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_dlstsq_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(dp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),                     intent(in)            :: b(:,:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp), allocatable, target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp), allocatable :: singular(:),rwork(:)
         real(dp), pointer :: xmat(:,:),amat(:,:)
         real(dp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call dgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_dlstsq_multiple


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_qlstsq_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(qp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),                     intent(in)            :: b(:,:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp), allocatable, target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp), allocatable :: singular(:),rwork(:)
         real(qp), pointer :: xmat(:,:),amat(:,:)
         real(qp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call qgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_qlstsq_multiple


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_clstsq_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),                     intent(in)            :: b(:,:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp), allocatable, target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp), allocatable :: singular(:),rwork(:)
         complex(sp), pointer :: xmat(:,:),amat(:,:)
         complex(sp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call cgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_clstsq_multiple


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_zlstsq_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),                     intent(in)            :: b(:,:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp), allocatable, target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp), allocatable :: singular(:),rwork(:)
         complex(dp), pointer :: xmat(:,:),amat(:,:)
         complex(dp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call zgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_zlstsq_multiple


     ! Compute the least-squares solution to a real system of linear equations Ax = B
     function stdlib_linalg_wlstsq_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),                     intent(inout), target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),                     intent(in)            :: b(:,:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp), allocatable, target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp), allocatable :: singular(:),rwork(:)
         complex(qp), pointer :: xmat(:,:),amat(:,:)
         complex(qp), allocatable :: cwork(:)
         character(*), parameter :: this = 'lstsq'

         !> Problem sizes
         m     = size(a,1,kind=ilp)
         lda   = size(a,1,kind=ilp)
         n     = size(a,2,kind=ilp)
         ldb   = size(b,1,kind=ilp)
         nrhs  = size(b  ,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda<1 .or. n<1 .or. ldb<1 .or. ldb/=m) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],',&
                                                                       'b=[',ldb,',',nrhs,']')
            allocate(x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(lda,n),source=a)
         else
            amat => a
         endif

         ! Initialize solution with the rhs
         allocate(x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in degreasing order)
         allocate(singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*S(1) are treated as zero.
         ! Use same default value as NumPy
         rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call wgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate(rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size a=[',lda,',',n,'], b[',ldb,',',nrhs,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not.copy_a) deallocate(amat)

         ! Process output and return
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_wlstsq_multiple


     ! Simple integer log2 implementation
     elemental integer(ilp) function ilog2(x)
        integer(ilp), intent(in) :: x

        integer(ilp) :: remndr

        if (x>0) then
           remndr = x
           ilog2 = -1_ilp
           do while (remndr>0)
               ilog2  = ilog2 + 1_ilp
               remndr = shiftr(remndr,1)
           end do
        else
           ilog2 = -huge(0_ilp)
        endif
     end function ilog2

end module stdlib_linalg_least_squares
