!!!#:include "common.fypp"
module stdlib_linalg_svd
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Singular values
     public :: svd

     !> LAPACK tasks

     !> Return full matrices U, V^T to separate storage
     character, parameter :: GESDD_FULL_MATRICES   = 'A'

     !> Return shrunk matrices U, V^T to k = min(m,n)
     character, parameter :: GESDD_SHRINK_MATRICES = 'S'

     !> Overwrite A storage with U (if M>=N) or VT (if M<N); separate storage for the other matrix
     character, parameter :: GESDD_OVERWRITE_A     = 'O'

     !> Do not return either U or VT (singular values array only)
     character, parameter :: GESDD_SINGVAL_ONLY    = 'N'

     character(*), parameter :: this = 'svd'

     ! Numpy: svd(a, full_matrices=True, compute_uv=True, hermitian=False)
     ! Scipy: svd(a, full_matrices=True, compute_uv=True, overwrite_a=False, check_finite=True, lapack_driver='gesdd')

     interface svd
        !!!#:for rk,rt,ri in ALL_KINDS_TYPES
        module procedure stdlib_linalg_svd_dp
        !!!#:endfor
     end interface svd


     contains

     !> Process GESDD output flag
     elemental subroutine gesdd_info(err,info,m,n)
        !> Error handler
        type(linalg_state), intent(inout) :: err
        !> GESDD return flag
        integer(ilp), intent(in) :: info
        !> Input matrix size
        integer(ilp), intent(in) :: m,n

        select case (info)
           case (0)
               ! Success!
               err%state = LINALG_SUCCESS
           case (-1)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid task ID on input to GESDD.')
           case (-5,-3:-2)
               err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',m,',',n,']')
           case (-8)
               err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix U size, with a=[',m,',',n,']')
           case (-10)
               err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix V size, with a=[',m,',',n,']')
           case (-4)
               err = linalg_state(this,LINALG_VALUE_ERROR,'A contains invalid/NaN values.')
           case (1:)
               err = linalg_state(this,LINALG_ERROR,'SVD computation did not converge.')
           case default
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Unknown error returned by GESDD.')
        end select

     end subroutine gesdd_info


     !!!#:for rk,rt,ri in ALL_KINDS_TYPES
     ! Return SVD of matrix A = U S V^T
     subroutine stdlib_linalg_svd_dp(a,s,u,vt,overwrite_a,full_matrices,err)
         !> Input matrix A[m,n]
         real(dp), intent(inout), target :: a(:,:)
         !> Array of singular values
         real(dp), intent(out) :: s(:)
         !> The columns of U contain the eigenvectors of A A^T
         real(dp), optional, intent(out), target :: u(:,:)
         !> The rows of V^T contain the eigenvectors fo A^T A
         real(dp), optional, intent(out), target :: vt(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk), optional, intent(in) :: overwrite_a
         !> [optional] full matrices have shape(u)==[m,m], shape(vh)==[n,n] (default); otherwise
         !> they are shape(u)==[m,k] and shape(vh)==[k,n] with k=min(m,n)
         logical(lk), optional, intent(in) :: full_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state), optional, intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldvt,info,k,lwork,liwork
         integer(ilp), allocatable :: iwork(:)
         logical(lk) :: copy_a,full_storage,compute_uv,alloc_u,alloc_vt,can_overwrite_a
         character :: task
         real(dp), target :: work_dummy(1),u_dummy(1,1),vt_dummy(1,1)
         real(dp), allocatable :: rwork(:)
         real(dp), pointer :: amat(:,:),umat(:,:),vtmat(:,:)

         !> Matrix determinant size
         m   = size(a,1,kind=ilp)
         n   = size(a,2,kind=ilp)
         k   = min(m,n)
         lda = m

         if (.not.k>0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size: a=[',m,',',n,']')
            goto 1
         end if

         if (.not.size(s,kind=ilp)>=k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'singular value array has insufficient size:',&
                                                        ' s=[',size(s,kind=ilp),'], k=',k)
            goto 1
         endif

         ! Integer storage
         liwork = 8*k
         allocate(iwork(liwork))

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not.overwrite_a
         else
            copy_a = .true._lk
         endif

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate(amat(m,n),source=a)

            ! Check if we can overwrite A with data that will be lost
            can_overwrite_a = merge(.not.present(u),.not.present(vt),m>=n)

         else
            amat => a

            can_overwrite_a = .false.

         endif

         ! Full-size matrices
         if (present(full_matrices)) then
            full_storage = full_matrices
         else
            full_storage = .true.
         endif

         ! Decide if U, VT matrices should be computed
         compute_uv = present(u) .or. present(vt)

         ! U, VT storage
         if (present(u)) then
            umat    => u
            alloc_u = .false.
         elseif ((copy_a .and. m>=n) .or. .not.compute_uv) then
            ! U not wanted, and A can be overwritten: do not allocate
            umat    => u_dummy
            alloc_u = .false.
         elseif (.not.full_storage) then
            allocate(umat(m,k))
            alloc_u = .true.
         else
            allocate(umat(m,m))
            alloc_u = .true.
         end if

         if (present(vt)) then
            vtmat    => vt
            alloc_vt = .false.
         elseif ((copy_a .and. m<n) .or. .not.compute_uv) then
            ! amat can be overwritten, VT not wanted: VT is returned upon A
            vtmat    => vt_dummy
            alloc_vt = .false.
         elseif (.not.full_storage) then
            allocate(vtmat(k,n))
            alloc_vt = .true.
         else
            allocate(vtmat(n,n))
            alloc_vt = .true.
         end if

         ldu  = size(umat ,1,kind=ilp)
         ldvt = size(vtmat,1,kind=ilp)

         ! Decide SVD task
         if (.not.compute_uv) then
            task = GESDD_SINGVAL_ONLY
         elseif (can_overwrite_a) then
            ! A is a copy: we can overwrite its storage
            task = GESDD_OVERWRITE_A
         elseif (.not.full_storage) then
            task = GESDD_SHRINK_MATRICES
         else
            task = GESDD_FULL_MATRICES
         end if

         ! Compute workspace
         lwork = -1_ilp
         call gesdd(task,m,n,amat,lda,s,umat,ldu,vtmat,ldvt,work_dummy,lwork,iwork,info)
         call gesdd_info(err0,info,m,n)

         ! Compute SVD
         if (info==0) then

            !> Prepare working storage
            lwork = nint(work_dummy(1), kind=ilp)
            allocate(rwork(lwork))

            !> Compute SVD
            call gesdd(task,m,n,amat,lda,s,umat,ldu,vtmat,ldvt,rwork,lwork,iwork,info)
            call gesdd_info(err0,info,m,n)

         endif

         ! Finalize storage and process output flag
         if (copy_a)    deallocate(amat)
         if (alloc_u)   deallocate(umat)
         if (alloc_vt)  deallocate(vtmat)
         1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_svd_dp


!!!     #:endfor

end module stdlib_linalg_svd
