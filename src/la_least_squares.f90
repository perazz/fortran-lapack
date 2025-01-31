!> Least squares solution interface
module la_least_squares
     use la_constants
     use la_blas
     use la_lapack
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> @brief Compute a least squares solution to system \f$ A \cdot x = b \f$,
     !!  i.e. such that the 2-norm \f$ \|b - A \cdot x\| \f$ is minimized.
     !!
     !! This function computes the least-squares solution to a real system of linear equations:
     !!
     !! \f$ A \cdot x = b \f$
     !!
     !! where A is an `n x n` matrix and b is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution `x` is returned as an allocatable array.
     !!
     !! @param[in,out] a The input matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] b The right-hand side vector (`n`) or matrix (`n x nrhs`).
     !! @param[in] cond (Optional) A cutoff for rank evaluation: singular values \f$ s(i) \f$ such that
     !!                \f$ s(i) \leq \text{cond} \cdot \max(s) \f$ are considered zero.
     !! @param[in] overwrite_a (Optional) If true, A and B may be overwritten and destroyed. Default is false.
     !! @param[out] rank (Optional) The rank of the matrix A.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     !! @return Solution matrix `x` of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK least-squares solvers such as `[*GELSS](@ref la_lapack::gelss)`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A and B may be lost.
     !!
     public :: lstsq

     interface lstsq
        module procedure la_slstsq_one
        module procedure la_dlstsq_one
        module procedure la_qlstsq_one
        module procedure la_clstsq_one
        module procedure la_zlstsq_one
        module procedure la_wlstsq_one
        module procedure la_slstsq_multiple
        module procedure la_dlstsq_multiple
        module procedure la_qlstsq_multiple
        module procedure la_clstsq_multiple
        module procedure la_zlstsq_multiple
        module procedure la_wlstsq_multiple
     end interface lstsq

     contains

     !> Workspace needed by real(sp) gesv
     subroutine sgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'sgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 12*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + mnmin*nrhs + (smlsiz + 1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine sgesv_space

     !> Workspace needed by real(dp) gesv
     subroutine dgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'dgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 12*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + mnmin*nrhs + (smlsiz + 1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine dgesv_space

     !> Workspace needed by real(qp) gesv
     subroutine qgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'qgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 12*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + mnmin*nrhs + (smlsiz + 1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine qgesv_space

     !> Workspace needed by complex(sp) gesv
     subroutine cgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'cgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 10*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*(1 + nrhs) + 2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine cgesv_space

     !> Workspace needed by complex(dp) gesv
     subroutine zgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'zgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 10*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*(1 + nrhs) + 2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine zgesv_space

     !> Workspace needed by complex(qp) gesv
     subroutine wgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'wgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 10*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*(1 + nrhs) + 2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine wgesv_space

     !> Compute the least-squares solution to a real(sp) system of linear equations Ax = B
     function la_slstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp),allocatable :: singular(:),rwork(:)
         real(sp),pointer :: xmat(:,:),amat(:,:)
         real(sp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call sgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_slstsq_one

     !> Compute the least-squares solution to a real(dp) system of linear equations Ax = B
     function la_dlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp),allocatable :: singular(:),rwork(:)
         real(dp),pointer :: xmat(:,:),amat(:,:)
         real(dp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call dgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_dlstsq_one

     !> Compute the least-squares solution to a real(qp) system of linear equations Ax = B
     function la_qlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp),allocatable :: singular(:),rwork(:)
         real(qp),pointer :: xmat(:,:),amat(:,:)
         real(qp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call qgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_qlstsq_one

     !> Compute the least-squares solution to a complex(sp) system of linear equations Ax = B
     function la_clstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp),allocatable :: singular(:),rwork(:)
         complex(sp),pointer :: xmat(:,:),amat(:,:)
         complex(sp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call cgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_clstsq_one

     !> Compute the least-squares solution to a complex(dp) system of linear equations Ax = B
     function la_zlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp),allocatable :: singular(:),rwork(:)
         complex(dp),pointer :: xmat(:,:),amat(:,:)
         complex(dp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call zgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_zlstsq_one

     !> Compute the least-squares solution to a complex(qp) system of linear equations Ax = B
     function la_wlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp),allocatable :: singular(:),rwork(:)
         complex(qp),pointer :: xmat(:,:),amat(:,:)
         complex(qp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call wgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_wlstsq_one

     !> Compute the least-squares solution to a real(sp) system of linear equations Ax = B
     function la_slstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp),allocatable :: singular(:),rwork(:)
         real(sp),pointer :: xmat(:,:),amat(:,:)
         real(sp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call sgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_slstsq_multiple

     !> Compute the least-squares solution to a real(dp) system of linear equations Ax = B
     function la_dlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp),allocatable :: singular(:),rwork(:)
         real(dp),pointer :: xmat(:,:),amat(:,:)
         real(dp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call dgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_dlstsq_multiple

     !> Compute the least-squares solution to a real(qp) system of linear equations Ax = B
     function la_qlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp),allocatable :: singular(:),rwork(:)
         real(qp),pointer :: xmat(:,:),amat(:,:)
         real(qp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call qgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,rwork,lrwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_qlstsq_multiple

     !> Compute the least-squares solution to a complex(sp) system of linear equations Ax = B
     function la_clstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(sp) :: acond,rcond
         real(sp),allocatable :: singular(:),rwork(:)
         complex(sp),pointer :: xmat(:,:),amat(:,:)
         complex(sp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Allocate working space
         call cgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_clstsq_multiple

     !> Compute the least-squares solution to a complex(dp) system of linear equations Ax = B
     function la_zlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(dp) :: acond,rcond
         real(dp),allocatable :: singular(:),rwork(:)
         complex(dp),pointer :: xmat(:,:),amat(:,:)
         complex(dp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Allocate working space
         call zgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_zlstsq_multiple

     !> Compute the least-squares solution to a complex(qp) system of linear equations Ax = B
     function la_wlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,nrhs,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp),allocatable :: iwork(:)
         logical(lk) :: copy_a
         real(qp) :: acond,rcond
         real(qp),allocatable :: singular(:),rwork(:)
         complex(qp),pointer :: xmat(:,:),amat(:,:)
         complex(qp),allocatable :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         mnmin = min(m,n)
         mnmax = max(m,n)

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                   'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            arank = 0
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Singular values array (in decreasing order)
         allocate (singular(mnmin))

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Allocate working space
         call wgesv_space(m,n,nrhs,lrwork,liwork,lcwork)
         allocate (rwork(lrwork),cwork(lcwork),iwork(liwork))

         ! Solve system using singular value decomposition
         call gelsd(m,n,nrhs,amat,lda,xmat,ldb,singular,rcond,arank,cwork,lrwork,rwork,iwork,info)

         ! The condition number of A in the 2-norm = S(1)/S(min(m,n)).
         acond = singular(1)/singular(mnmin)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
            case (1:)
                err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call err0%handle(err)
         if (present(rank)) rank = arank

     end function la_wlstsq_multiple

     ! Simple integer log2 implementation
     elemental integer(ilp) function ilog2(x)
        integer(ilp),intent(in) :: x

        integer(ilp) :: remndr

        if (x > 0) then
           remndr = x
           ilog2 = -1_ilp
           do while (remndr > 0)
               ilog2 = ilog2 + 1_ilp
               remndr = shiftr(remndr,1)
           end do
        else
           ilog2 = -huge(0_ilp)
        end if
     end function ilog2

end module la_least_squares
