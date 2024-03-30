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

     ! Numpy: eigenvalues, eigenvectors = eig(a)
     !        eigenvalues = eigvals(a)
     ! Scipy: eig(a, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True, homogeneous_eigvals=False)

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
        module procedure stdlib_linalg_eigvals_d
        module procedure stdlib_linalg_eigvals_q
        module procedure stdlib_linalg_eigvals_c
        module procedure stdlib_linalg_eigvals_z
        module procedure stdlib_linalg_eigvals_w
     end interface eigvals

     character(*),parameter :: this = 'eigenvalues'

     contains
     
     !> Request for eigenvector calculation
     elemental character function eigenvectors_flag(required)
        logical,intent(in) :: required
        eigenvectors_flag = merge('V','N',required)
     end function eigenvectors_flag

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
               err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',m,',',n,']')
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

     !> Singular values of matrix A
     function stdlib_linalg_eigvals_s(a,err) result(lambda)
         !> Input matrix A[m,n]
         real(sp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
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

     !> SVD of matrix A = U S V^T, returning S and optionally U and V
     subroutine stdlib_linalg_eig_s(a,lambda,left,right,overwrite_a,err)
         !> Input matrix A[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(sp),intent(out) :: lambda(:)
         !> The columns of LEFT contain the left eigenvectors of A
         real(sp),optional,intent(out),target :: left(:,:)
         !> The columns of RIGHT contain the right eigenvectors of A
         real(sp),optional,intent(out),target :: right(:,:)
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
         real(sp),pointer :: amat(:,:),umat(:,:),vmat(:,:),lreal(:),limag(:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n],', must be square.')
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
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace
         allocate (lreal(n),limag(n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lreal,limag, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,info)
         call geev_info(err0,info,m,n)

         ! Compute SVD
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
         
2        if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_s

     !> Singular values of matrix A
     function stdlib_linalg_eigvals_d(a,err) result(lambda)
         !> Input matrix A[m,n]
         real(dp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
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

     !> SVD of matrix A = U S V^T, returning S and optionally U and V
     subroutine stdlib_linalg_eig_d(a,lambda,left,right,overwrite_a,err)
         !> Input matrix A[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(dp),intent(out) :: lambda(:)
         !> The columns of LEFT contain the left eigenvectors of A
         real(dp),optional,intent(out),target :: left(:,:)
         !> The columns of RIGHT contain the right eigenvectors of A
         real(dp),optional,intent(out),target :: right(:,:)
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
         real(dp),pointer :: amat(:,:),umat(:,:),vmat(:,:),lreal(:),limag(:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n],', must be square.')
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
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace
         allocate (lreal(n),limag(n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lreal,limag, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,info)
         call geev_info(err0,info,m,n)

         ! Compute SVD
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
         
2        if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_d

     !> Singular values of matrix A
     function stdlib_linalg_eigvals_q(a,err) result(lambda)
         !> Input matrix A[m,n]
         real(qp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
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

     !> SVD of matrix A = U S V^T, returning S and optionally U and V
     subroutine stdlib_linalg_eig_q(a,lambda,left,right,overwrite_a,err)
         !> Input matrix A[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(qp),intent(out) :: lambda(:)
         !> The columns of LEFT contain the left eigenvectors of A
         real(qp),optional,intent(out),target :: left(:,:)
         !> The columns of RIGHT contain the right eigenvectors of A
         real(qp),optional,intent(out),target :: right(:,:)
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
         real(qp),pointer :: amat(:,:),umat(:,:),vmat(:,:),lreal(:),limag(:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n],', must be square.')
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
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace
         allocate (lreal(n),limag(n))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lreal,limag, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,info)
         call geev_info(err0,info,m,n)

         ! Compute SVD
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
         
2        if (copy_a) deallocate (amat)
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_eig_q

     !> Singular values of matrix A
     function stdlib_linalg_eigvals_c(a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(sp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
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

     !> SVD of matrix A = U S V^T, returning S and optionally U and V
     subroutine stdlib_linalg_eig_c(a,lambda,left,right,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(sp),intent(out) :: lambda(:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(sp),optional,intent(out),target :: left(:,:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(sp),optional,intent(out),target :: right(:,:)
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
         complex(sp),pointer :: amat(:,:),umat(:,:),vmat(:,:),lreal(:),limag(:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n],', must be square.')
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
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace
!         if (task==geev_SINGVAL_ONLY) then
!            lrwork = max(1,7*k)
!         else
!            lrwork = max(1,5*k*(k+1),2*k*(k+max(m,n))+k)
!         endif
!         allocate(rwork(lrwork))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lambda, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,rwork,info)
         call geev_info(err0,info,m,n)

         ! Compute SVD
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

     !> Singular values of matrix A
     function stdlib_linalg_eigvals_z(a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(dp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
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

     !> SVD of matrix A = U S V^T, returning S and optionally U and V
     subroutine stdlib_linalg_eig_z(a,lambda,left,right,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(dp),intent(out) :: lambda(:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(dp),optional,intent(out),target :: left(:,:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(dp),optional,intent(out),target :: right(:,:)
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
         complex(dp),pointer :: amat(:,:),umat(:,:),vmat(:,:),lreal(:),limag(:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n],', must be square.')
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
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace
!         if (task==geev_SINGVAL_ONLY) then
!            lrwork = max(1,7*k)
!         else
!            lrwork = max(1,5*k*(k+1),2*k*(k+max(m,n))+k)
!         endif
!         allocate(rwork(lrwork))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lambda, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,rwork,info)
         call geev_info(err0,info,m,n)

         ! Compute SVD
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

     !> Singular values of matrix A
     function stdlib_linalg_eigvals_w(a,err) result(lambda)
         !> Input matrix A[m,n]
         complex(qp),intent(in),target :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
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

     !> SVD of matrix A = U S V^T, returning S and optionally U and V
     subroutine stdlib_linalg_eig_w(a,lambda,left,right,overwrite_a,err)
         !> Input matrix A[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Array of eigenvalues
         complex(qp),intent(out) :: lambda(:)
         !> The columns of LEFT contain the left eigenvectors of A
         complex(qp),optional,intent(out),target :: left(:,:)
         !> The columns of RIGHT contain the right eigenvectors of A
         complex(qp),optional,intent(out),target :: right(:,:)
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
         complex(qp),pointer :: amat(:,:),umat(:,:),vmat(:,:),lreal(:),limag(:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         neig = size(lambda,kind=ilp)
         lda = m

         if (.not. (k > 0 .and. m == n)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or matrix size a=', [m,n],', must be square.')
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
                                        shape(vmat),', with n=',n)
               goto 2
            end if
            
         else
            umat => u_dummy
         end if

         ldu = size(umat,1,kind=ilp)
         ldv = size(vmat,1,kind=ilp)

         ! Compute workspace
!         if (task==geev_SINGVAL_ONLY) then
!            lrwork = max(1,7*k)
!         else
!            lrwork = max(1,5*k*(k+1),2*k*(k+max(m,n))+k)
!         endif
!         allocate(rwork(lrwork))

         lwork = -1_ilp
        
         call geev(task_u,task_v,n,amat,lda, &
                   lambda, &
                   umat,ldu,vmat,ldv, &
                   work_dummy,lwork,rwork,info)
         call geev_info(err0,info,m,n)

         ! Compute SVD
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

end module stdlib_linalg_eig
