module stdlib_linalg_eig
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Eigendecomposition of a square matrix: return eigenvalues, and optionally eigenvectors
     public :: eig
     !> Eigenvalues of a square matrix
     public :: eigvals
     !> Eigendecomposition of a real symmetric or complex hermitian matrix
     public :: eigh
     !> Eigenvalues of a real symmetric or complex hermitian matrix
     public :: eigvalsh

     ! Numpy: eigenvalues, eigenvectors = eig(a)
     !        eigenvalues = eigvals(a)
     ! Scipy: eig(a, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True, homogeneous_eigvals=False)

     ! Numpy: eigenvalues, eigenvectors = eigh(a, uplo='L')
     !        eigenvalues = eigvalsh(a)
     ! Scipy: eigh(a, b=None, *, lower=True, eigvals_only=False, overwrite_a=False, overwrite_b=False, turbo=<object object>, eigvals=<object object>, type=1, check_finite=True, subset_by_index=None, subset_by_value=None, driver=None)

     interface eig
        module procedure stdlib_linalg_eig_s
        module procedure stdlib_linalg_eig_d
        module procedure stdlib_linalg_eig_q
        module procedure stdlib_linalg_eig_c
        module procedure stdlib_linalg_eig_z
        module procedure stdlib_linalg_eig_w
     end interface eig

     interface eigvals
        module procedure stdlib_linalg_eigvals_s
        module procedure stdlib_linalg_eigvals_noerr_s
        module procedure stdlib_linalg_eigvals_d
        module procedure stdlib_linalg_eigvals_noerr_d
        module procedure stdlib_linalg_eigvals_q
        module procedure stdlib_linalg_eigvals_noerr_q
        module procedure stdlib_linalg_eigvals_c
        module procedure stdlib_linalg_eigvals_noerr_c
        module procedure stdlib_linalg_eigvals_z
        module procedure stdlib_linalg_eigvals_noerr_z
        module procedure stdlib_linalg_eigvals_w
        module procedure stdlib_linalg_eigvals_noerr_w
     end interface eigvals
     
     interface eigh
        module procedure stdlib_linalg_eigh_s
        module procedure stdlib_linalg_eigh_d
        module procedure stdlib_linalg_eigh_q
        module procedure stdlib_linalg_eigh_c
        module procedure stdlib_linalg_eigh_z
        module procedure stdlib_linalg_eigh_w
     end interface eigh
     
     interface eigvalsh
        module procedure stdlib_linalg_eigvalsh_s
        module procedure stdlib_linalg_eigvalsh_noerr_s
        module procedure stdlib_linalg_eigvalsh_d
        module procedure stdlib_linalg_eigvalsh_noerr_d
        module procedure stdlib_linalg_eigvalsh_q
        module procedure stdlib_linalg_eigvalsh_noerr_q
        module procedure stdlib_linalg_eigvalsh_c
        module procedure stdlib_linalg_eigvalsh_noerr_c
        module procedure stdlib_linalg_eigvalsh_z
        module procedure stdlib_linalg_eigvalsh_noerr_z
        module procedure stdlib_linalg_eigvalsh_w
        module procedure stdlib_linalg_eigvalsh_noerr_w
     end interface eigvalsh

     character(*),parameter :: this = 'eigenvalues'

     contains
     
     !> Request for eigenvector calculation
     elemental character function eigenvectors_flag(required)
        logical(lk),intent(in) :: required
        eigenvectors_flag = merge('V','N',required)
     end function eigenvectors_flag
     
     !> Request for symmetry side (default: lower)
     elemental character function symmetric_triangle(upper)
        logical(lk),optional,intent(in) :: upper
        if (present(upper)) then
           symmetric_triangle = merge('U','L',upper)
        else
           symmetric_triangle = 'L'
        end if
     end function symmetric_triangle

     !> Process GEEV output flags
     elemental subroutine geev_info(err,info,m,n)
        !> Error handler
        type(linalg_state),intent(inout) :: err
        !> geev return flag
        integer(ilp),intent(in) :: info
        !> Input matrix size
        integer(ilp),intent(in) :: m,n

        select case (info)
           case (0)
               ! Success!
               err%state = LINALG_SUCCESS
           case (-1)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid task ID: left eigenvectors.')
           case (-2)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid task ID: right eigenvectors.')
           case (-5,-3)
               err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
           case (-9)
               err = linalg_state(this,LINALG_VALUE_ERROR,'insufficient left vector matrix size.')
           case (-11)
               err = linalg_state(this,LINALG_VALUE_ERROR,'insufficient right vector matrix size.')
           case (-13)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Insufficient work array size.')
           case (1:)
               err = linalg_state(this,LINALG_ERROR,'Eigenvalue computation did not converge.')
           case default
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Unknown error returned by geev.')
        end select

     end subroutine geev_info

     !> Process SYEV/HEEV output flags
     elemental subroutine heev_info(err,info,m,n)
        !> Error handler
        type(linalg_state),intent(inout) :: err
        !> geev return flag
        integer(ilp),intent(in) :: info
        !> Input matrix size
        integer(ilp),intent(in) :: m,n

        select case (info)
           case (0)
               ! Success!
               err%state = LINALG_SUCCESS
           case (-1)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid eigenvector request.')
           case (-2)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid triangular section request.')
           case (-5,-3)
               err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
           case (-8)
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'insufficient workspace size.')
           case (1:)
               err = linalg_state(this,LINALG_ERROR,'Eigenvalue computation did not converge.')
           case default
               err = linalg_state(this,LINALG_INTERNAL_ERROR,'Unknown error returned by syev/heev.')
        end select

     end subroutine heev_info

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_s(a,err) result(lambda)
         !> Input matrix A[m,n]
         real(sp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         complex(sp),allocatable :: lambda(:)

         !> Create
         real(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_s(amat,lambda,overwrite_a=.false.,err=err)

     end function stdlib_linalg_eigvals_s

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_noerr_s(a) result(lambda)
         !> Input matrix A[m,n]
         real(sp),intent(in),target :: a(:,:)
         !> Array of singular values
         complex(sp),allocatable :: lambda(:)

         !> Create
         real(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_s(amat,lambda,overwrite_a=.false.)

     end function stdlib_linalg_eigvals_noerr_s

     !> Eigendecomposition of matrix A returning an array `lambda` of eigenvalues,
     !> and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eig_s(a,lambda,right,left,overwrite_a,err)
         !> Input matrix A[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(sp),intent(out) :: lambda(:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(sp),optional,intent(out),target :: right(:,:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(sp),optional,intent(out),target :: left(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldv,info,k,lwork,lrwork,neig
         logical(lk) :: copy_a
         character :: task_u,task_v
         real(sp),target :: work_dummy(1),u_dummy(1,1),v_dummy(1,1)
         real(sp),allocatable :: work(:)
         real(sp),allocatable :: rwork(:)
         real(sp),pointer :: amat(:,:),lreal(:),limag(:),umat(:,:),vmat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,', n=',n)
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
            allocate (amat(m,n),source=a)
         else
            amat => a
         end if

         ! Decide if U, V eigenvectors
         task_u = eigenvectors_flag(present(left))
         task_v = eigenvectors_flag(present(right))
         
         if (present(right)) then
            
            allocate (vmat(n,n))
            
            if (size(vmat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'right eigenvector matrix has insufficient size: ', &
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            vmat => v_dummy
         end if
            
         if (present(left)) then
            
            allocate (umat(n,n))
            
            if (size(umat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'left eigenvector matrix has insufficient size: ', &
                                        shape(umat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace size
         allocate (lreal(n),limag(n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lreal,limag, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,info)
         call geev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=sp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call geev(task_u,task_v,n,amat,lda, &
                      lreal,limag, &
                      umat,ldu,vmat,ldv, &
                      work,lwork,info)
            call geev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         lambda(:n) = cmplx(lreal(:n),limag(:n),kind=sp)
         
         ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
         ! geev returns reals as:
         ! u(j)   = VL(:,j) + i*VL(:,j+1) and
         ! u(j+1) = VL(:,j) - i*VL(:,j+1).
         ! Convert these to complex numbers here.
         if (present(right)) call assign_real_eigenvectors_sp(n,lambda,vmat,right)
         if (present(left)) call assign_real_eigenvectors_sp(n,lambda,umat,left)
         
2        if (copy_a) deallocate (amat)
         if (present(right)) deallocate (vmat)
         if (present(left)) deallocate (umat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_s
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_s(a,upper_a,err) result(lambda)
         !> Input matrix A[m,n]
         real(sp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         real(sp),allocatable :: lambda(:)
         
         real(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_s(amat,lambda,upper_a=upper_a,overwrite_a=.false.,err=err)
         
     end function stdlib_linalg_eigvalsh_s
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_noerr_s(a,upper_a) result(lambda)
         !> Input matrix A[m,n]
         real(sp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> Array of singular values
         real(sp),allocatable :: lambda(:)

         real(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_s(amat,lambda,upper_a=upper_a,overwrite_a=.false.)

     end function stdlib_linalg_eigvalsh_noerr_s

     !> Eigendecomposition of a real symmetric or complex Hermitian matrix A returning an array `lambda`
     !> of eigenvalues, and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eigh_s(a,lambda,vectors,upper_a,overwrite_a,err)
         !> Input matrix A[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         real(sp),intent(out) :: lambda(:)
         !> The columns of vectors contain the orthonormal eigenvectors of A
         real(sp),optional,intent(out),target :: vectors(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,info,k,lwork,neig
         logical(lk) :: copy_a
         character :: triangle,task
         real(sp),target :: work_dummy(1)
         real(sp),allocatable :: work(:)
         real(sp),allocatable :: rwork(:)
         real(sp),pointer :: amat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,' must be >= n=',n)
            goto 1
         end if
        
         ! Check if input A can be overwritten
         if (present(vectors)) then
            ! No need to copy A anyways
            copy_a = .false.
         elseif (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if
         
         ! Should we use the upper or lower half of the matrix?
         triangle = symmetric_triangle(upper_a)
         
         ! Request for eigenvectors
         task = eigenvectors_flag(present(vectors))
         
         if (present(vectors)) then
            
            ! Check size
            if (any(shape(vectors,kind=ilp) < n)) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'eigenvector matrix has insufficient size: ', &
                                        shape(vectors),', with n=',n)
               goto 1
            end if
            
            ! The input matrix will be overwritten by the vectors.
            ! So, use this one as storage for syev/heev
            amat => vectors
            
            ! Copy data in
            amat(:n,:n) = a(:n,:n)
                        
         elseif (copy_a) then
            ! Initialize a matrix temporary
            allocate (amat(m,n),source=a)
         else
            ! Overwrite A
            amat => a
         end if

         lda = size(amat,1,kind=ilp)

         ! Request workspace size
         lwork = -1_ilp
         call syev(task,triangle,n,amat,lda,lambda,work_dummy,lwork,info)
         call heev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=sp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call syev(task,triangle,n,amat,lda,lambda,work,lwork,info)
            call heev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eigh_s
     
     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_d(a,err) result(lambda)
         !> Input matrix A[m,n]
         real(dp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         complex(dp),allocatable :: lambda(:)

         !> Create
         real(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_d(amat,lambda,overwrite_a=.false.,err=err)

     end function stdlib_linalg_eigvals_d

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_noerr_d(a) result(lambda)
         !> Input matrix A[m,n]
         real(dp),intent(in),target :: a(:,:)
         !> Array of singular values
         complex(dp),allocatable :: lambda(:)

         !> Create
         real(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_d(amat,lambda,overwrite_a=.false.)

     end function stdlib_linalg_eigvals_noerr_d

     !> Eigendecomposition of matrix A returning an array `lambda` of eigenvalues,
     !> and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eig_d(a,lambda,right,left,overwrite_a,err)
         !> Input matrix A[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(dp),intent(out) :: lambda(:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(dp),optional,intent(out),target :: right(:,:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(dp),optional,intent(out),target :: left(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldv,info,k,lwork,lrwork,neig
         logical(lk) :: copy_a
         character :: task_u,task_v
         real(dp),target :: work_dummy(1),u_dummy(1,1),v_dummy(1,1)
         real(dp),allocatable :: work(:)
         real(dp),allocatable :: rwork(:)
         real(dp),pointer :: amat(:,:),lreal(:),limag(:),umat(:,:),vmat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,', n=',n)
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
            allocate (amat(m,n),source=a)
         else
            amat => a
         end if

         ! Decide if U, V eigenvectors
         task_u = eigenvectors_flag(present(left))
         task_v = eigenvectors_flag(present(right))
         
         if (present(right)) then
            
            allocate (vmat(n,n))
            
            if (size(vmat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'right eigenvector matrix has insufficient size: ', &
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            vmat => v_dummy
         end if
            
         if (present(left)) then
            
            allocate (umat(n,n))
            
            if (size(umat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'left eigenvector matrix has insufficient size: ', &
                                        shape(umat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace size
         allocate (lreal(n),limag(n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lreal,limag, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,info)
         call geev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=dp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call geev(task_u,task_v,n,amat,lda, &
                      lreal,limag, &
                      umat,ldu,vmat,ldv, &
                      work,lwork,info)
            call geev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         lambda(:n) = cmplx(lreal(:n),limag(:n),kind=dp)
         
         ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
         ! geev returns reals as:
         ! u(j)   = VL(:,j) + i*VL(:,j+1) and
         ! u(j+1) = VL(:,j) - i*VL(:,j+1).
         ! Convert these to complex numbers here.
         if (present(right)) call assign_real_eigenvectors_dp(n,lambda,vmat,right)
         if (present(left)) call assign_real_eigenvectors_dp(n,lambda,umat,left)
         
2        if (copy_a) deallocate (amat)
         if (present(right)) deallocate (vmat)
         if (present(left)) deallocate (umat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_d
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_d(a,upper_a,err) result(lambda)
         !> Input matrix A[m,n]
         real(dp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         real(dp),allocatable :: lambda(:)
         
         real(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_d(amat,lambda,upper_a=upper_a,overwrite_a=.false.,err=err)
         
     end function stdlib_linalg_eigvalsh_d
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_noerr_d(a,upper_a) result(lambda)
         !> Input matrix A[m,n]
         real(dp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> Array of singular values
         real(dp),allocatable :: lambda(:)

         real(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_d(amat,lambda,upper_a=upper_a,overwrite_a=.false.)

     end function stdlib_linalg_eigvalsh_noerr_d

     !> Eigendecomposition of a real symmetric or complex Hermitian matrix A returning an array `lambda`
     !> of eigenvalues, and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eigh_d(a,lambda,vectors,upper_a,overwrite_a,err)
         !> Input matrix A[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         real(dp),intent(out) :: lambda(:)
         !> The columns of vectors contain the orthonormal eigenvectors of A
         real(dp),optional,intent(out),target :: vectors(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,info,k,lwork,neig
         logical(lk) :: copy_a
         character :: triangle,task
         real(dp),target :: work_dummy(1)
         real(dp),allocatable :: work(:)
         real(dp),allocatable :: rwork(:)
         real(dp),pointer :: amat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,' must be >= n=',n)
            goto 1
         end if
        
         ! Check if input A can be overwritten
         if (present(vectors)) then
            ! No need to copy A anyways
            copy_a = .false.
         elseif (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if
         
         ! Should we use the upper or lower half of the matrix?
         triangle = symmetric_triangle(upper_a)
         
         ! Request for eigenvectors
         task = eigenvectors_flag(present(vectors))
         
         if (present(vectors)) then
            
            ! Check size
            if (any(shape(vectors,kind=ilp) < n)) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'eigenvector matrix has insufficient size: ', &
                                        shape(vectors),', with n=',n)
               goto 1
            end if
            
            ! The input matrix will be overwritten by the vectors.
            ! So, use this one as storage for syev/heev
            amat => vectors
            
            ! Copy data in
            amat(:n,:n) = a(:n,:n)
                        
         elseif (copy_a) then
            ! Initialize a matrix temporary
            allocate (amat(m,n),source=a)
         else
            ! Overwrite A
            amat => a
         end if

         lda = size(amat,1,kind=ilp)

         ! Request workspace size
         lwork = -1_ilp
         call syev(task,triangle,n,amat,lda,lambda,work_dummy,lwork,info)
         call heev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=dp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call syev(task,triangle,n,amat,lda,lambda,work,lwork,info)
            call heev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eigh_d
     
     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_q(a,err) result(lambda)
         !> Input matrix A[m,n]
         real(qp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         complex(qp),allocatable :: lambda(:)

         !> Create
         real(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_q(amat,lambda,overwrite_a=.false.,err=err)

     end function stdlib_linalg_eigvals_q

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_noerr_q(a) result(lambda)
         !> Input matrix A[m,n]
         real(qp),intent(in),target :: a(:,:)
         !> Array of singular values
         complex(qp),allocatable :: lambda(:)

         !> Create
         real(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_q(amat,lambda,overwrite_a=.false.)

     end function stdlib_linalg_eigvals_noerr_q

     !> Eigendecomposition of matrix A returning an array `lambda` of eigenvalues,
     !> and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eig_q(a,lambda,right,left,overwrite_a,err)
         !> Input matrix A[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(qp),intent(out) :: lambda(:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(qp),optional,intent(out),target :: right(:,:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(qp),optional,intent(out),target :: left(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldv,info,k,lwork,lrwork,neig
         logical(lk) :: copy_a
         character :: task_u,task_v
         real(qp),target :: work_dummy(1),u_dummy(1,1),v_dummy(1,1)
         real(qp),allocatable :: work(:)
         real(qp),allocatable :: rwork(:)
         real(qp),pointer :: amat(:,:),lreal(:),limag(:),umat(:,:),vmat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,', n=',n)
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
            allocate (amat(m,n),source=a)
         else
            amat => a
         end if

         ! Decide if U, V eigenvectors
         task_u = eigenvectors_flag(present(left))
         task_v = eigenvectors_flag(present(right))
         
         if (present(right)) then
            
            allocate (vmat(n,n))
            
            if (size(vmat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'right eigenvector matrix has insufficient size: ', &
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            vmat => v_dummy
         end if
            
         if (present(left)) then
            
            allocate (umat(n,n))
            
            if (size(umat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'left eigenvector matrix has insufficient size: ', &
                                        shape(umat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace size
         allocate (lreal(n),limag(n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lreal,limag, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,info)
         call geev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=qp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call geev(task_u,task_v,n,amat,lda, &
                      lreal,limag, &
                      umat,ldu,vmat,ldv, &
                      work,lwork,info)
            call geev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         lambda(:n) = cmplx(lreal(:n),limag(:n),kind=qp)
         
         ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
         ! geev returns reals as:
         ! u(j)   = VL(:,j) + i*VL(:,j+1) and
         ! u(j+1) = VL(:,j) - i*VL(:,j+1).
         ! Convert these to complex numbers here.
         if (present(right)) call assign_real_eigenvectors_qp(n,lambda,vmat,right)
         if (present(left)) call assign_real_eigenvectors_qp(n,lambda,umat,left)
         
2        if (copy_a) deallocate (amat)
         if (present(right)) deallocate (vmat)
         if (present(left)) deallocate (umat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_q
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_q(a,upper_a,err) result(lambda)
         !> Input matrix A[m,n]
         real(qp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         real(qp),allocatable :: lambda(:)
         
         real(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_q(amat,lambda,upper_a=upper_a,overwrite_a=.false.,err=err)
         
     end function stdlib_linalg_eigvalsh_q
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_noerr_q(a,upper_a) result(lambda)
         !> Input matrix A[m,n]
         real(qp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> Array of singular values
         real(qp),allocatable :: lambda(:)

         real(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_q(amat,lambda,upper_a=upper_a,overwrite_a=.false.)

     end function stdlib_linalg_eigvalsh_noerr_q

     !> Eigendecomposition of a real symmetric or complex Hermitian matrix A returning an array `lambda`
     !> of eigenvalues, and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eigh_q(a,lambda,vectors,upper_a,overwrite_a,err)
         !> Input matrix A[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         real(qp),intent(out) :: lambda(:)
         !> The columns of vectors contain the orthonormal eigenvectors of A
         real(qp),optional,intent(out),target :: vectors(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,info,k,lwork,neig
         logical(lk) :: copy_a
         character :: triangle,task
         real(qp),target :: work_dummy(1)
         real(qp),allocatable :: work(:)
         real(qp),allocatable :: rwork(:)
         real(qp),pointer :: amat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,' must be >= n=',n)
            goto 1
         end if
        
         ! Check if input A can be overwritten
         if (present(vectors)) then
            ! No need to copy A anyways
            copy_a = .false.
         elseif (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if
         
         ! Should we use the upper or lower half of the matrix?
         triangle = symmetric_triangle(upper_a)
         
         ! Request for eigenvectors
         task = eigenvectors_flag(present(vectors))
         
         if (present(vectors)) then
            
            ! Check size
            if (any(shape(vectors,kind=ilp) < n)) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'eigenvector matrix has insufficient size: ', &
                                        shape(vectors),', with n=',n)
               goto 1
            end if
            
            ! The input matrix will be overwritten by the vectors.
            ! So, use this one as storage for syev/heev
            amat => vectors
            
            ! Copy data in
            amat(:n,:n) = a(:n,:n)
                        
         elseif (copy_a) then
            ! Initialize a matrix temporary
            allocate (amat(m,n),source=a)
         else
            ! Overwrite A
            amat => a
         end if

         lda = size(amat,1,kind=ilp)

         ! Request workspace size
         lwork = -1_ilp
         call syev(task,triangle,n,amat,lda,lambda,work_dummy,lwork,info)
         call heev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=qp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call syev(task,triangle,n,amat,lda,lambda,work,lwork,info)
            call heev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eigh_q
     
     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_c(a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(sp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         complex(sp),allocatable :: lambda(:)

         !> Create
         complex(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_c(amat,lambda,overwrite_a=.false.,err=err)

     end function stdlib_linalg_eigvals_c

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_noerr_c(a) result(lambda)
         !> Input matrix A[m,n]
         complex(sp),intent(in),target :: a(:,:)
         !> Array of singular values
         complex(sp),allocatable :: lambda(:)

         !> Create
         complex(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_c(amat,lambda,overwrite_a=.false.)

     end function stdlib_linalg_eigvals_noerr_c

     !> Eigendecomposition of matrix A returning an array `lambda` of eigenvalues,
     !> and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eig_c(a,lambda,right,left,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(sp),intent(out) :: lambda(:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(sp),optional,intent(out),target :: right(:,:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(sp),optional,intent(out),target :: left(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldv,info,k,lwork,lrwork,neig
         logical(lk) :: copy_a
         character :: task_u,task_v
         complex(sp),target :: work_dummy(1),u_dummy(1,1),v_dummy(1,1)
         complex(sp),allocatable :: work(:)
         real(sp),allocatable :: rwork(:)
         complex(sp),pointer :: amat(:,:),lreal(:),limag(:),umat(:,:),vmat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,', n=',n)
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
            allocate (amat(m,n),source=a)
         else
            amat => a
         end if

         ! Decide if U, V eigenvectors
         task_u = eigenvectors_flag(present(left))
         task_v = eigenvectors_flag(present(right))
         
         if (present(right)) then
            
            vmat => right
            
            if (size(vmat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'right eigenvector matrix has insufficient size: ', &
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            vmat => v_dummy
         end if
            
         if (present(left)) then
            
            umat => left
            
            if (size(umat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'left eigenvector matrix has insufficient size: ', &
                                        shape(umat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace size
         allocate (rwork(2*n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lambda, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,rwork,info)
         call geev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=sp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call geev(task_u,task_v,n,amat,lda, &
                      lambda, &
                      umat,ldu,vmat,ldv, &
                      work,lwork,rwork,info)
            call geev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         
2        if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_c
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_c(a,upper_a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(sp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         real(sp),allocatable :: lambda(:)
         
         complex(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_c(amat,lambda,upper_a=upper_a,overwrite_a=.false.,err=err)
         
     end function stdlib_linalg_eigvalsh_c
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_noerr_c(a,upper_a) result(lambda)
         !> Input matrix A[m,n]
         complex(sp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> Array of singular values
         real(sp),allocatable :: lambda(:)

         complex(sp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_c(amat,lambda,upper_a=upper_a,overwrite_a=.false.)

     end function stdlib_linalg_eigvalsh_noerr_c

     !> Eigendecomposition of a real symmetric or complex Hermitian matrix A returning an array `lambda`
     !> of eigenvalues, and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eigh_c(a,lambda,vectors,upper_a,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         real(sp),intent(out) :: lambda(:)
         !> The columns of vectors contain the orthonormal eigenvectors of A
         complex(sp),optional,intent(out),target :: vectors(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,info,k,lwork,neig
         logical(lk) :: copy_a
         character :: triangle,task
         complex(sp),target :: work_dummy(1)
         complex(sp),allocatable :: work(:)
         real(sp),allocatable :: rwork(:)
         complex(sp),pointer :: amat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,' must be >= n=',n)
            goto 1
         end if
        
         ! Check if input A can be overwritten
         if (present(vectors)) then
            ! No need to copy A anyways
            copy_a = .false.
         elseif (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if
         
         ! Should we use the upper or lower half of the matrix?
         triangle = symmetric_triangle(upper_a)
         
         ! Request for eigenvectors
         task = eigenvectors_flag(present(vectors))
         
         if (present(vectors)) then
            
            ! Check size
            if (any(shape(vectors,kind=ilp) < n)) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'eigenvector matrix has insufficient size: ', &
                                        shape(vectors),', with n=',n)
               goto 1
            end if
            
            ! The input matrix will be overwritten by the vectors.
            ! So, use this one as storage for syev/heev
            amat => vectors
            
            ! Copy data in
            amat(:n,:n) = a(:n,:n)
                        
         elseif (copy_a) then
            ! Initialize a matrix temporary
            allocate (amat(m,n),source=a)
         else
            ! Overwrite A
            amat => a
         end if

         lda = size(amat,1,kind=ilp)

         ! Request workspace size
         lwork = -1_ilp
         allocate (rwork(max(1,3*n - 2)))
         call heev(task,triangle,n,amat,lda,lambda,work_dummy,lwork,rwork,info)
         call heev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=sp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call heev(task,triangle,n,amat,lda,lambda,work,lwork,rwork,info)
            call heev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eigh_c
     
     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_z(a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(dp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         complex(dp),allocatable :: lambda(:)

         !> Create
         complex(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_z(amat,lambda,overwrite_a=.false.,err=err)

     end function stdlib_linalg_eigvals_z

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_noerr_z(a) result(lambda)
         !> Input matrix A[m,n]
         complex(dp),intent(in),target :: a(:,:)
         !> Array of singular values
         complex(dp),allocatable :: lambda(:)

         !> Create
         complex(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_z(amat,lambda,overwrite_a=.false.)

     end function stdlib_linalg_eigvals_noerr_z

     !> Eigendecomposition of matrix A returning an array `lambda` of eigenvalues,
     !> and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eig_z(a,lambda,right,left,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(dp),intent(out) :: lambda(:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(dp),optional,intent(out),target :: right(:,:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(dp),optional,intent(out),target :: left(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldv,info,k,lwork,lrwork,neig
         logical(lk) :: copy_a
         character :: task_u,task_v
         complex(dp),target :: work_dummy(1),u_dummy(1,1),v_dummy(1,1)
         complex(dp),allocatable :: work(:)
         real(dp),allocatable :: rwork(:)
         complex(dp),pointer :: amat(:,:),lreal(:),limag(:),umat(:,:),vmat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,', n=',n)
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
            allocate (amat(m,n),source=a)
         else
            amat => a
         end if

         ! Decide if U, V eigenvectors
         task_u = eigenvectors_flag(present(left))
         task_v = eigenvectors_flag(present(right))
         
         if (present(right)) then
            
            vmat => right
            
            if (size(vmat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'right eigenvector matrix has insufficient size: ', &
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            vmat => v_dummy
         end if
            
         if (present(left)) then
            
            umat => left
            
            if (size(umat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'left eigenvector matrix has insufficient size: ', &
                                        shape(umat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace size
         allocate (rwork(2*n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lambda, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,rwork,info)
         call geev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=dp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call geev(task_u,task_v,n,amat,lda, &
                      lambda, &
                      umat,ldu,vmat,ldv, &
                      work,lwork,rwork,info)
            call geev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         
2        if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_z
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_z(a,upper_a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(dp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         real(dp),allocatable :: lambda(:)
         
         complex(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_z(amat,lambda,upper_a=upper_a,overwrite_a=.false.,err=err)
         
     end function stdlib_linalg_eigvalsh_z
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_noerr_z(a,upper_a) result(lambda)
         !> Input matrix A[m,n]
         complex(dp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> Array of singular values
         real(dp),allocatable :: lambda(:)

         complex(dp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_z(amat,lambda,upper_a=upper_a,overwrite_a=.false.)

     end function stdlib_linalg_eigvalsh_noerr_z

     !> Eigendecomposition of a real symmetric or complex Hermitian matrix A returning an array `lambda`
     !> of eigenvalues, and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eigh_z(a,lambda,vectors,upper_a,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         real(dp),intent(out) :: lambda(:)
         !> The columns of vectors contain the orthonormal eigenvectors of A
         complex(dp),optional,intent(out),target :: vectors(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,info,k,lwork,neig
         logical(lk) :: copy_a
         character :: triangle,task
         complex(dp),target :: work_dummy(1)
         complex(dp),allocatable :: work(:)
         real(dp),allocatable :: rwork(:)
         complex(dp),pointer :: amat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,' must be >= n=',n)
            goto 1
         end if
        
         ! Check if input A can be overwritten
         if (present(vectors)) then
            ! No need to copy A anyways
            copy_a = .false.
         elseif (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if
         
         ! Should we use the upper or lower half of the matrix?
         triangle = symmetric_triangle(upper_a)
         
         ! Request for eigenvectors
         task = eigenvectors_flag(present(vectors))
         
         if (present(vectors)) then
            
            ! Check size
            if (any(shape(vectors,kind=ilp) < n)) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'eigenvector matrix has insufficient size: ', &
                                        shape(vectors),', with n=',n)
               goto 1
            end if
            
            ! The input matrix will be overwritten by the vectors.
            ! So, use this one as storage for syev/heev
            amat => vectors
            
            ! Copy data in
            amat(:n,:n) = a(:n,:n)
                        
         elseif (copy_a) then
            ! Initialize a matrix temporary
            allocate (amat(m,n),source=a)
         else
            ! Overwrite A
            amat => a
         end if

         lda = size(amat,1,kind=ilp)

         ! Request workspace size
         lwork = -1_ilp
         allocate (rwork(max(1,3*n - 2)))
         call heev(task,triangle,n,amat,lda,lambda,work_dummy,lwork,rwork,info)
         call heev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=dp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call heev(task,triangle,n,amat,lda,lambda,work,lwork,rwork,info)
            call heev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eigh_z
     
     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_w(a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(qp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         complex(qp),allocatable :: lambda(:)

         !> Create
         complex(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_w(amat,lambda,overwrite_a=.false.,err=err)

     end function stdlib_linalg_eigvals_w

     !> Return an array of eigenvalues of matrix A.
     function stdlib_linalg_eigvals_noerr_w(a) result(lambda)
         !> Input matrix A[m,n]
         complex(qp),intent(in),target :: a(:,:)
         !> Array of singular values
         complex(qp),allocatable :: lambda(:)

         !> Create
         complex(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eig_w(amat,lambda,overwrite_a=.false.)

     end function stdlib_linalg_eigvals_noerr_w

     !> Eigendecomposition of matrix A returning an array `lambda` of eigenvalues,
     !> and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eig_w(a,lambda,right,left,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(qp),intent(out) :: lambda(:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(qp),optional,intent(out),target :: right(:,:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(qp),optional,intent(out),target :: left(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,ldu,ldv,info,k,lwork,lrwork,neig
         logical(lk) :: copy_a
         character :: task_u,task_v
         complex(qp),target :: work_dummy(1),u_dummy(1,1),v_dummy(1,1)
         complex(qp),allocatable :: work(:)
         real(qp),allocatable :: rwork(:)
         complex(qp),pointer :: amat(:,:),lreal(:),limag(:),umat(:,:),vmat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,', n=',n)
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
            allocate (amat(m,n),source=a)
         else
            amat => a
         end if

         ! Decide if U, V eigenvectors
         task_u = eigenvectors_flag(present(left))
         task_v = eigenvectors_flag(present(right))
         
         if (present(right)) then
            
            vmat => right
            
            if (size(vmat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'right eigenvector matrix has insufficient size: ', &
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            vmat => v_dummy
         end if
            
         if (present(left)) then
            
            umat => left
            
            if (size(umat,2,kind=ilp) < n) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'left eigenvector matrix has insufficient size: ', &
                                        shape(umat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace size
         allocate (rwork(2*n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lambda, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,rwork,info)
         call geev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=qp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call geev(task_u,task_v,n,amat,lda, &
                      lambda, &
                      umat,ldu,vmat,ldv, &
                      work,lwork,rwork,info)
            call geev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         
2        if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_w
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_w(a,upper_a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(qp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Array of singular values
         real(qp),allocatable :: lambda(:)
         
         complex(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_w(amat,lambda,upper_a=upper_a,overwrite_a=.false.,err=err)
         
     end function stdlib_linalg_eigvalsh_w
     
     !> Return an array of eigenvalues of real symmetric / complex hermitian A
     function stdlib_linalg_eigvalsh_noerr_w(a,upper_a) result(lambda)
         !> Input matrix A[m,n]
         complex(qp),intent(in),target :: a(:,:)
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> Array of singular values
         real(qp),allocatable :: lambda(:)

         complex(qp),pointer :: amat(:,:)
         integer(ilp) :: m,n,k

         !> Create an internal pointer so the intent of A won't affect the next call
         amat => a

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)

         !> Allocate return storage
         allocate (lambda(k))

         !> Compute eigenvalues only
         call stdlib_linalg_eigh_w(amat,lambda,upper_a=upper_a,overwrite_a=.false.)

     end function stdlib_linalg_eigvalsh_noerr_w

     !> Eigendecomposition of a real symmetric or complex Hermitian matrix A returning an array `lambda`
     !> of eigenvalues, and optionally right or left eigenvectors.
     subroutine stdlib_linalg_eigh_w(a,lambda,vectors,upper_a,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         real(qp),intent(out) :: lambda(:)
         !> The columns of vectors contain the orthonormal eigenvectors of A
         complex(qp),optional,intent(out),target :: vectors(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Should the upper/lower half of A be used? Default: lower
         logical(lk),optional,intent(in) :: upper_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,lda,info,k,lwork,neig
         logical(lk) :: copy_a
         character :: triangle,task
         complex(qp),target :: work_dummy(1)
         complex(qp),allocatable :: work(:)
         real(qp),allocatable :: rwork(:)
         complex(qp),pointer :: amat(:,:)

         !> Matrix size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n], &
                                                        ', must be square.')
            goto 1
         end if

         if (.not. neig >= k) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'eigenvalue array has insufficient size:', &
                                                        ' lambda=',neig,' must be >= n=',n)
            goto 1
         end if
        
         ! Check if input A can be overwritten
         if (present(vectors)) then
            ! No need to copy A anyways
            copy_a = .false.
         elseif (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if
         
         ! Should we use the upper or lower half of the matrix?
         triangle = symmetric_triangle(upper_a)
         
         ! Request for eigenvectors
         task = eigenvectors_flag(present(vectors))
         
         if (present(vectors)) then
            
            ! Check size
            if (any(shape(vectors,kind=ilp) < n)) then
               err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                                        'eigenvector matrix has insufficient size: ', &
                                        shape(vectors),', with n=',n)
               goto 1
            end if
            
            ! The input matrix will be overwritten by the vectors.
            ! So, use this one as storage for syev/heev
            amat => vectors
            
            ! Copy data in
            amat(:n,:n) = a(:n,:n)
                        
         elseif (copy_a) then
            ! Initialize a matrix temporary
            allocate (amat(m,n),source=a)
         else
            ! Overwrite A
            amat => a
         end if

         lda = size(amat,1,kind=ilp)

         ! Request workspace size
         lwork = -1_ilp
         allocate (rwork(max(1,3*n - 2)))
         call heev(task,triangle,n,amat,lda,lambda,work_dummy,lwork,rwork,info)
         call heev_info(err0,info,m,n)

         ! Compute eigenvalues
         if (info == 0) then

            !> Prepare working storage
            lwork = nint(real(work_dummy(1),kind=qp),kind=ilp)
            allocate (work(lwork))

            !> Compute eigensystem
            call heev(task,triangle,n,amat,lda,lambda,work,lwork,rwork,info)
            call heev_info(err0,info,m,n)

         end if
         
         ! Finalize storage and process output flag
         if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eigh_w
     
     !> GEEV for real matrices returns complex eigenvalues in real arrays.
     !> Convert them to complex here, following the GEEV logic.
     pure subroutine assign_real_eigenvectors_sp(n,lambda,lmat,out_mat)
        !> Problem size
        integer(ilp),intent(in) :: n
        !> Array of eigenvalues
        complex(sp),intent(in) :: lambda(:)
        !> Real matrix as returned by geev
        real(sp),intent(in) :: lmat(:,:)
        !> Complex matrix as returned by eig
        complex(sp),intent(out) :: out_mat(:,:)
        
        integer(ilp) :: i,j
        
        ! Copy matrix
        do concurrent(i=1:n,j=1:n)
           out_mat(i,j) = lmat(i,j)
        end do
        
        ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
        ! geev returns them as reals as:
        ! u(j)   = VL(:,j) + i*VL(:,j+1) and
        ! u(j+1) = VL(:,j) - i*VL(:,j+1).
        ! Convert these to complex numbers here.
        do j = 1,n - 1
           if (lambda(j) == conjg(lambda(j + 1))) then
              out_mat(:,j) = cmplx(lmat(:,j),lmat(:,j + 1),kind=sp)
              out_mat(:,j + 1) = cmplx(lmat(:,j),-lmat(:,j + 1),kind=sp)
           end if
        end do
        
     end subroutine assign_real_eigenvectors_sp
     !> GEEV for real matrices returns complex eigenvalues in real arrays.
     !> Convert them to complex here, following the GEEV logic.
     pure subroutine assign_real_eigenvectors_dp(n,lambda,lmat,out_mat)
        !> Problem size
        integer(ilp),intent(in) :: n
        !> Array of eigenvalues
        complex(dp),intent(in) :: lambda(:)
        !> Real matrix as returned by geev
        real(dp),intent(in) :: lmat(:,:)
        !> Complex matrix as returned by eig
        complex(dp),intent(out) :: out_mat(:,:)
        
        integer(ilp) :: i,j
        
        ! Copy matrix
        do concurrent(i=1:n,j=1:n)
           out_mat(i,j) = lmat(i,j)
        end do
        
        ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
        ! geev returns them as reals as:
        ! u(j)   = VL(:,j) + i*VL(:,j+1) and
        ! u(j+1) = VL(:,j) - i*VL(:,j+1).
        ! Convert these to complex numbers here.
        do j = 1,n - 1
           if (lambda(j) == conjg(lambda(j + 1))) then
              out_mat(:,j) = cmplx(lmat(:,j),lmat(:,j + 1),kind=dp)
              out_mat(:,j + 1) = cmplx(lmat(:,j),-lmat(:,j + 1),kind=dp)
           end if
        end do
        
     end subroutine assign_real_eigenvectors_dp
     !> GEEV for real matrices returns complex eigenvalues in real arrays.
     !> Convert them to complex here, following the GEEV logic.
     pure subroutine assign_real_eigenvectors_qp(n,lambda,lmat,out_mat)
        !> Problem size
        integer(ilp),intent(in) :: n
        !> Array of eigenvalues
        complex(qp),intent(in) :: lambda(:)
        !> Real matrix as returned by geev
        real(qp),intent(in) :: lmat(:,:)
        !> Complex matrix as returned by eig
        complex(qp),intent(out) :: out_mat(:,:)
        
        integer(ilp) :: i,j
        
        ! Copy matrix
        do concurrent(i=1:n,j=1:n)
           out_mat(i,j) = lmat(i,j)
        end do
        
        ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
        ! geev returns them as reals as:
        ! u(j)   = VL(:,j) + i*VL(:,j+1) and
        ! u(j+1) = VL(:,j) - i*VL(:,j+1).
        ! Convert these to complex numbers here.
        do j = 1,n - 1
           if (lambda(j) == conjg(lambda(j + 1))) then
              out_mat(:,j) = cmplx(lmat(:,j),lmat(:,j + 1),kind=qp)
              out_mat(:,j + 1) = cmplx(lmat(:,j),-lmat(:,j + 1),kind=qp)
           end if
        end do
        
     end subroutine assign_real_eigenvectors_qp

end module stdlib_linalg_eig
