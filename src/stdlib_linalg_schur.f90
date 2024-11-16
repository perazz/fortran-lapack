module stdlib_linalg_schur
    use stdlib_linalg_constants
    use stdlib_linalg_lapack,only:gees
    use stdlib_linalg_state,only:linalg_state,linalg_error_handling,LINALG_ERROR, &
        LINALG_INTERNAL_ERROR,LINALG_VALUE_ERROR
    implicit none(type,external)
    private
    
    public :: schur
    public :: schur_space

    character(*),parameter :: this = 'schur'
    
    !> List of internal GEES tasks:
    
    !> No task request
    character,parameter :: GEES_NOT = 'N'
    
    !> Request Schur vectors to be computed
    character,parameter :: GEES_WITH_VECTORS = 'V'
    
    !> Request Schur vectors to be sorted
    character,parameter :: GEES_SORTED_VECTORS = 'S'
    
  ! Schur decomposition of rank-2 array A
  interface schur
    !! version: experimental
    !!
    !! Computes the Schur decomposition of matrix \( A = Z T Z^H \).
    !! ([Specification](../page/specs/stdlib_linalg.html#schur-compute-the-schur-decomposition-of-a-matrix))
    !!
    !!### Summary
    !! Compute the Schur decomposition of a `real` or `complex` matrix: \( A = Z T Z^H \), where \( Z \) is
    !! orthonormal/unitary and \( T \) is upper-triangular or quasi-upper-triangular. Matrix \( A \) has size `[m,m]`.
    !!
    !!### Description
    !!
    !! This interface provides methods for computing the Schur decomposition of a matrix.
    !! Supported data types include `real` and `complex`. If a pre-allocated workspace is provided, no internal
    !! memory allocations take place when using this interface.
    !!
    !! The output matrix \( T \) is upper-triangular for `complex` input matrices and quasi-upper-triangular
    !! for `real` input matrices (with possible \( 2 \times 2 \) blocks on the diagonal).
    !! The user can optionally request sorting of eigenvalues based on conditions, or a custom sorting function.
    !!
    !!@note The solution is based on LAPACK's Schur decomposition routines (`*GEES`). Sorting options
    !! are implemented using LAPACK's eigenvalue sorting mechanism.
    !!
      module procedure stdlib_linalg_s_schur
      module procedure stdlib_linalg_d_schur
      module procedure stdlib_linalg_q_schur
      module procedure stdlib_linalg_c_schur
      module procedure stdlib_linalg_z_schur
      module procedure stdlib_linalg_w_schur
  end interface schur

  ! Return the working array space required by the Schur decomposition solver
  interface schur_space
    !! version: experimental
    !!
    !! Computes the working array space required by the Schur decomposition solver
    !! ([Specification](../page/specs/stdlib_linalg.html#schur-space-compute-internal-working-space-requirements-for-the-schur-decomposition))
    !!
    !!### Description
    !!
    !! This interface returns the size of the `real` or `complex` working storage required by the
    !! Schur decomposition solver. The working size only depends on the kind (`real` or `complex`) and size of
    !! the matrix being decomposed. Storage size can be used to pre-allocate a working array in case several
    !! repeated Schur decompositions of same-size matrices are sought. If pre-allocated working arrays
    !! are provided, no internal allocations will take place during the decomposition.
    !!
      module procedure get_schur_s_workspace
      module procedure get_schur_d_workspace
      module procedure get_schur_q_workspace
      module procedure get_schur_c_workspace
      module procedure get_schur_z_workspace
      module procedure get_schur_w_workspace
  end interface schur_space
    
    contains

    !> Wrapper function for Schur vectors request
    elemental character function gees_vectors(wanted)
        !> Are Schur vectors wanted?
        logical(lk),intent(in) :: wanted
        gees_vectors = merge(GEES_WITH_VECTORS,GEES_NOT,wanted)
    end function gees_vectors
    
    !> Wrapper function for Schur vectors request
    elemental character function gees_sort_eigs(sorted)
        !> Should the eigenvalues be sorted?
        logical(lk),intent(in) :: sorted
        gees_sort_eigs = merge(GEES_SORTED_VECTORS,GEES_NOT,sorted)
    end function gees_sort_eigs
    
    !> Wrapper function to handle GEES error codes
    elemental subroutine handle_gees_info(info,m,n,ldvs,err)
        integer(ilp),intent(in) :: info,m,n,ldvs
        type(linalg_state),intent(out) :: err

        ! Process GEES output
        select case (info)
        case (0_ilp)
            ! Success
        case (-1_ilp)
            ! Vector not wanted, but task is wrong
            err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid Schur vector task request')
        case (-2_ilp)
            ! Vector not wanted, but task is wrong
            err = linalg_state(this,LINALG_INTERNAL_ERROR,'Invalid sorting task request')
        case (-4_ilp,-6_ilp)
            ! Vector not wanted, but task is wrong
            err = linalg_state(this,LINALG_VALUE_ERROR,'Invalid/non-square input matrix size:', [m,n])
        case (-11_ilp)
            err = linalg_state(this,LINALG_VALUE_ERROR,'Schur vector matrix has insufficient size', [ldvs,n])
        case (-13_ilp)
            err = linalg_state(this,LINALG_INTERNAL_ERROR,'Insufficient working storage size')
        case (1_ilp:)
            
            if (info == n + 2) then
                err = linalg_state(this,LINALG_ERROR,'Ill-conditioned problem: could not sort eigenvalues')
            elseif (info == n + 1) then
                err = linalg_state(this,LINALG_ERROR,'Some selected eigenvalues lost property due to sorting')
            else
                err = linalg_state(this,LINALG_ERROR,'Convergence failure; converged range is', [info,n])
            end if
            
        case default
            
            err = linalg_state(this,LINALG_INTERNAL_ERROR,'GEES catastrophic error: info=',info)

        end select
        
    end subroutine handle_gees_info
    
    subroutine get_schur_s_workspace(a,lwork,err)
        !> Input matrix a[m,m]
        real(sp),intent(in),target :: a(:,:)
        !> Minimum workspace size for the decomposition operation
        integer(ilp),intent(out) :: lwork
        !> State return flag. Returns an error if the query failed
        type(linalg_state),optional,intent(out) :: err
        
        integer(ilp) :: m,n,sdim,info
        character :: jobvs,sort
        logical(lk) :: bwork_dummy(1)
        real(sp),pointer :: amat(:,:)
        real(sp) :: rwork_dummy(1)
        real(sp) :: wr_dummy(1),wi_dummy(1),vs_dummy(1,1),work_dummy(1)
        type(linalg_state) :: err0
        
        !> Initialize problem
        lwork = -1_ilp
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        
        !> Create a dummy intent(inout) argument
        amat => a
        
        !> Select dummy task
        jobvs = gees_vectors(.true.)
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        !> Get Schur workspace
        call gees(jobvs,sort,do_not_select,n,amat,m,sdim,wr_dummy,wi_dummy, &
                  vs_dummy,m,work_dummy,lwork,bwork_dummy,info)
        if (info == 0) lwork = int(work_dummy(1),kind=ilp)
        call handle_gees_info(info,m,n,m,err0)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Interface to dummy select routine
            pure logical(lk) function do_not_select(alphar,alphai)
                real(sp),intent(in) :: alphar,alphai
                do_not_select = .false.
            end function do_not_select
        
    end subroutine get_schur_s_workspace
    
    ! Schur decomposition subroutine
    pure subroutine stdlib_linalg_s_schur(a,t,z,storage,err)
        !> Input matrix a[m,m]
        real(sp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(sp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(sp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(sp),intent(out),optional,target :: storage(:)
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
!        type(linalg_state) :: err0
!        integer(ilp) :: m, lda, info, liwork
!        logical(lk) :: overwrite_a_
!        logical, pointer :: bwork(:)
!        integer(ilp), allocatable :: iwork(:)
!        real(sp), pointer :: amat(:,:), tau(:), work(:)
!        real(sp) :: rwork_dummy(1)  ! Dummy for real/complex cases
!        real(sp), allocatable :: tmat(:,:), zmat(:,:)
!        character :: jobz

!        ! Problem size
!        m = size(a, 1, kind=ilp)
!        lda = size(a, 1, kind=ilp)
!
!        ! Validate dimensions
!        if (size(a, 1, kind=ilp) /= size(a, 2, kind=ilp)) then
!            err0 = linalg_state(this, LINALG_VALUE_ERROR, 'Matrix A must be square: a=', [size(a,1), size(a,2)])
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Set default values
!        overwrite_a_ = .false._lk
!        if (present(overwrite_a)) overwrite_a_ = overwrite_a
!
!        ! Job type based on sorting request
!        if (present(sort) .and. sort) then
!            jobz = 'S'  ! Compute and sort eigenvalues
!        else
!            jobz = 'N'  ! Compute Schur decomposition without sorting
!        end if
!
!        ! Prepare storage
!        allocate(tmat(m, m), source=0.0_sp)
!        allocate(zmat(m, m), source=0.0_sp)
!
!        if (overwrite_a_) then
!            amat => a
!        else
!            allocate(amat(m, m), source=a)
!        end if
!
!        ! Allocate workspace
!        liwork = -1
!        if (present(lwork)) then
!            allocate(work(lwork))
!        else
!            allocate(work(1))  ! Temporary workspace for querying size
!        end if
!
!        ! Workspace query
!        call  gees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, .false._lk, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Update workspace size
!        if (.not.present(lwork)) then
!            liwork = ceiling(real(work(1), kind=sp), kind=ilp)
!            deallocate(work)
!            allocate(work(liwork))
!        end if
!
!        ! Compute Schur decomposition
!        call  gees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, present(sort) .and. sort, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Output results
!        t = amat
!        z = zmat

    end subroutine stdlib_linalg_s_schur

    subroutine get_schur_d_workspace(a,lwork,err)
        !> Input matrix a[m,m]
        real(dp),intent(in),target :: a(:,:)
        !> Minimum workspace size for the decomposition operation
        integer(ilp),intent(out) :: lwork
        !> State return flag. Returns an error if the query failed
        type(linalg_state),optional,intent(out) :: err
        
        integer(ilp) :: m,n,sdim,info
        character :: jobvs,sort
        logical(lk) :: bwork_dummy(1)
        real(dp),pointer :: amat(:,:)
        real(dp) :: rwork_dummy(1)
        real(dp) :: wr_dummy(1),wi_dummy(1),vs_dummy(1,1),work_dummy(1)
        type(linalg_state) :: err0
        
        !> Initialize problem
        lwork = -1_ilp
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        
        !> Create a dummy intent(inout) argument
        amat => a
        
        !> Select dummy task
        jobvs = gees_vectors(.true.)
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        !> Get Schur workspace
        call gees(jobvs,sort,do_not_select,n,amat,m,sdim,wr_dummy,wi_dummy, &
                  vs_dummy,m,work_dummy,lwork,bwork_dummy,info)
        if (info == 0) lwork = int(work_dummy(1),kind=ilp)
        call handle_gees_info(info,m,n,m,err0)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Interface to dummy select routine
            pure logical(lk) function do_not_select(alphar,alphai)
                real(dp),intent(in) :: alphar,alphai
                do_not_select = .false.
            end function do_not_select
        
    end subroutine get_schur_d_workspace
    
    ! Schur decomposition subroutine
    pure subroutine stdlib_linalg_d_schur(a,t,z,storage,err)
        !> Input matrix a[m,m]
        real(dp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(dp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(dp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(dp),intent(out),optional,target :: storage(:)
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
!        type(linalg_state) :: err0
!        integer(ilp) :: m, lda, info, liwork
!        logical(lk) :: overwrite_a_
!        logical, pointer :: bwork(:)
!        integer(ilp), allocatable :: iwork(:)
!        real(dp), pointer :: amat(:,:), tau(:), work(:)
!        real(dp) :: rwork_dummy(1)  ! Dummy for real/complex cases
!        real(dp), allocatable :: tmat(:,:), zmat(:,:)
!        character :: jobz

!        ! Problem size
!        m = size(a, 1, kind=ilp)
!        lda = size(a, 1, kind=ilp)
!
!        ! Validate dimensions
!        if (size(a, 1, kind=ilp) /= size(a, 2, kind=ilp)) then
!            err0 = linalg_state(this, LINALG_VALUE_ERROR, 'Matrix A must be square: a=', [size(a,1), size(a,2)])
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Set default values
!        overwrite_a_ = .false._lk
!        if (present(overwrite_a)) overwrite_a_ = overwrite_a
!
!        ! Job type based on sorting request
!        if (present(sort) .and. sort) then
!            jobz = 'S'  ! Compute and sort eigenvalues
!        else
!            jobz = 'N'  ! Compute Schur decomposition without sorting
!        end if
!
!        ! Prepare storage
!        allocate(tmat(m, m), source=0.0_dp)
!        allocate(zmat(m, m), source=0.0_dp)
!
!        if (overwrite_a_) then
!            amat => a
!        else
!            allocate(amat(m, m), source=a)
!        end if
!
!        ! Allocate workspace
!        liwork = -1
!        if (present(lwork)) then
!            allocate(work(lwork))
!        else
!            allocate(work(1))  ! Temporary workspace for querying size
!        end if
!
!        ! Workspace query
!        call  gees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, .false._lk, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Update workspace size
!        if (.not.present(lwork)) then
!            liwork = ceiling(real(work(1), kind=dp), kind=ilp)
!            deallocate(work)
!            allocate(work(liwork))
!        end if
!
!        ! Compute Schur decomposition
!        call  gees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, present(sort) .and. sort, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Output results
!        t = amat
!        z = zmat

    end subroutine stdlib_linalg_d_schur

    subroutine get_schur_q_workspace(a,lwork,err)
        !> Input matrix a[m,m]
        real(qp),intent(in),target :: a(:,:)
        !> Minimum workspace size for the decomposition operation
        integer(ilp),intent(out) :: lwork
        !> State return flag. Returns an error if the query failed
        type(linalg_state),optional,intent(out) :: err
        
        integer(ilp) :: m,n,sdim,info
        character :: jobvs,sort
        logical(lk) :: bwork_dummy(1)
        real(qp),pointer :: amat(:,:)
        real(qp) :: rwork_dummy(1)
        real(qp) :: wr_dummy(1),wi_dummy(1),vs_dummy(1,1),work_dummy(1)
        type(linalg_state) :: err0
        
        !> Initialize problem
        lwork = -1_ilp
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        
        !> Create a dummy intent(inout) argument
        amat => a
        
        !> Select dummy task
        jobvs = gees_vectors(.true.)
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        !> Get Schur workspace
        call gees(jobvs,sort,do_not_select,n,amat,m,sdim,wr_dummy,wi_dummy, &
                  vs_dummy,m,work_dummy,lwork,bwork_dummy,info)
        if (info == 0) lwork = int(work_dummy(1),kind=ilp)
        call handle_gees_info(info,m,n,m,err0)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Interface to dummy select routine
            pure logical(lk) function do_not_select(alphar,alphai)
                real(qp),intent(in) :: alphar,alphai
                do_not_select = .false.
            end function do_not_select
        
    end subroutine get_schur_q_workspace
    
    ! Schur decomposition subroutine
    pure subroutine stdlib_linalg_q_schur(a,t,z,storage,err)
        !> Input matrix a[m,m]
        real(qp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(qp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(qp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(qp),intent(out),optional,target :: storage(:)
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
!        type(linalg_state) :: err0
!        integer(ilp) :: m, lda, info, liwork
!        logical(lk) :: overwrite_a_
!        logical, pointer :: bwork(:)
!        integer(ilp), allocatable :: iwork(:)
!        real(qp), pointer :: amat(:,:), tau(:), work(:)
!        real(qp) :: rwork_dummy(1)  ! Dummy for real/complex cases
!        real(qp), allocatable :: tmat(:,:), zmat(:,:)
!        character :: jobz

!        ! Problem size
!        m = size(a, 1, kind=ilp)
!        lda = size(a, 1, kind=ilp)
!
!        ! Validate dimensions
!        if (size(a, 1, kind=ilp) /= size(a, 2, kind=ilp)) then
!            err0 = linalg_state(this, LINALG_VALUE_ERROR, 'Matrix A must be square: a=', [size(a,1), size(a,2)])
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Set default values
!        overwrite_a_ = .false._lk
!        if (present(overwrite_a)) overwrite_a_ = overwrite_a
!
!        ! Job type based on sorting request
!        if (present(sort) .and. sort) then
!            jobz = 'S'  ! Compute and sort eigenvalues
!        else
!            jobz = 'N'  ! Compute Schur decomposition without sorting
!        end if
!
!        ! Prepare storage
!        allocate(tmat(m, m), source=0.0_qp)
!        allocate(zmat(m, m), source=0.0_qp)
!
!        if (overwrite_a_) then
!            amat => a
!        else
!            allocate(amat(m, m), source=a)
!        end if
!
!        ! Allocate workspace
!        liwork = -1
!        if (present(lwork)) then
!            allocate(work(lwork))
!        else
!            allocate(work(1))  ! Temporary workspace for querying size
!        end if
!
!        ! Workspace query
!        call  gees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, .false._lk, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Update workspace size
!        if (.not.present(lwork)) then
!            liwork = ceiling(real(work(1), kind=qp), kind=ilp)
!            deallocate(work)
!            allocate(work(liwork))
!        end if
!
!        ! Compute Schur decomposition
!        call  gees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, present(sort) .and. sort, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Output results
!        t = amat
!        z = zmat

    end subroutine stdlib_linalg_q_schur

    subroutine get_schur_c_workspace(a,lwork,err)
        !> Input matrix a[m,m]
        complex(sp),intent(in),target :: a(:,:)
        !> Minimum workspace size for the decomposition operation
        integer(ilp),intent(out) :: lwork
        !> State return flag. Returns an error if the query failed
        type(linalg_state),optional,intent(out) :: err
        
        integer(ilp) :: m,n,sdim,info
        character :: jobvs,sort
        logical(lk) :: bwork_dummy(1)
        complex(sp),pointer :: amat(:,:)
        real(sp) :: rwork_dummy(1)
        complex(sp) :: wr_dummy(1),wi_dummy(1),vs_dummy(1,1),work_dummy(1)
        type(linalg_state) :: err0
        
        !> Initialize problem
        lwork = -1_ilp
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        
        !> Create a dummy intent(inout) argument
        amat => a
        
        !> Select dummy task
        jobvs = gees_vectors(.true.)
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        !> Get Schur workspace
        call gees(jobvs,sort,do_not_select,n,amat,m,sdim,wr_dummy, &
                  vs_dummy,m,work_dummy,lwork,rwork_dummy,bwork_dummy,info)
        if (info == 0) lwork = int(work_dummy(1),kind=ilp)
        call handle_gees_info(info,m,n,m,err0)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Interface to dummy select routine
            pure logical(lk) function do_not_select(alpha)
                complex(sp),intent(in) :: alpha
                do_not_select = .false.
            end function do_not_select
        
    end subroutine get_schur_c_workspace
    
    ! Schur decomposition subroutine
    pure subroutine stdlib_linalg_c_schur(a,t,z,storage,err)
        !> Input matrix a[m,m]
        complex(sp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(sp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(sp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(sp),intent(out),optional,target :: storage(:)
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
!        type(linalg_state) :: err0
!        integer(ilp) :: m, lda, info, liwork
!        logical(lk) :: overwrite_a_
!        logical, pointer :: bwork(:)
!        integer(ilp), allocatable :: iwork(:)
!        complex(sp), pointer :: amat(:,:), tau(:), work(:)
!        complex(sp) :: rwork_dummy(1)  ! Dummy for real/complex cases
!        complex(sp), allocatable :: tmat(:,:), zmat(:,:)
!        character :: jobz

!        ! Problem size
!        m = size(a, 1, kind=ilp)
!        lda = size(a, 1, kind=ilp)
!
!        ! Validate dimensions
!        if (size(a, 1, kind=ilp) /= size(a, 2, kind=ilp)) then
!            err0 = linalg_state(this, LINALG_VALUE_ERROR, 'Matrix A must be square: a=', [size(a,1), size(a,2)])
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Set default values
!        overwrite_a_ = .false._lk
!        if (present(overwrite_a)) overwrite_a_ = overwrite_a
!
!        ! Job type based on sorting request
!        if (present(sort) .and. sort) then
!            jobz = 'S'  ! Compute and sort eigenvalues
!        else
!            jobz = 'N'  ! Compute Schur decomposition without sorting
!        end if
!
!        ! Prepare storage
!        allocate(tmat(m, m), source=0.0_sp)
!        allocate(zmat(m, m), source=0.0_sp)
!
!        if (overwrite_a_) then
!            amat => a
!        else
!            allocate(amat(m, m), source=a)
!        end if
!
!        ! Allocate workspace
!        liwork = -1
!        if (present(lwork)) then
!            allocate(work(lwork))
!        else
!            allocate(work(1))  ! Temporary workspace for querying size
!        end if
!
!        ! Workspace query
!        call  zgees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, .false._lk, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Update workspace size
!        if (.not.present(lwork)) then
!            liwork = ceiling(real(work(1), kind=sp), kind=ilp)
!            deallocate(work)
!            allocate(work(liwork))
!        end if
!
!        ! Compute Schur decomposition
!        call  zgees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, present(sort) .and. sort, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Output results
!        t = amat
!        z = zmat

    end subroutine stdlib_linalg_c_schur

    subroutine get_schur_z_workspace(a,lwork,err)
        !> Input matrix a[m,m]
        complex(dp),intent(in),target :: a(:,:)
        !> Minimum workspace size for the decomposition operation
        integer(ilp),intent(out) :: lwork
        !> State return flag. Returns an error if the query failed
        type(linalg_state),optional,intent(out) :: err
        
        integer(ilp) :: m,n,sdim,info
        character :: jobvs,sort
        logical(lk) :: bwork_dummy(1)
        complex(dp),pointer :: amat(:,:)
        real(dp) :: rwork_dummy(1)
        complex(dp) :: wr_dummy(1),wi_dummy(1),vs_dummy(1,1),work_dummy(1)
        type(linalg_state) :: err0
        
        !> Initialize problem
        lwork = -1_ilp
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        
        !> Create a dummy intent(inout) argument
        amat => a
        
        !> Select dummy task
        jobvs = gees_vectors(.true.)
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        !> Get Schur workspace
        call gees(jobvs,sort,do_not_select,n,amat,m,sdim,wr_dummy, &
                  vs_dummy,m,work_dummy,lwork,rwork_dummy,bwork_dummy,info)
        if (info == 0) lwork = int(work_dummy(1),kind=ilp)
        call handle_gees_info(info,m,n,m,err0)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Interface to dummy select routine
            pure logical(lk) function do_not_select(alpha)
                complex(dp),intent(in) :: alpha
                do_not_select = .false.
            end function do_not_select
        
    end subroutine get_schur_z_workspace
    
    ! Schur decomposition subroutine
    pure subroutine stdlib_linalg_z_schur(a,t,z,storage,err)
        !> Input matrix a[m,m]
        complex(dp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(dp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(dp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(dp),intent(out),optional,target :: storage(:)
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
!        type(linalg_state) :: err0
!        integer(ilp) :: m, lda, info, liwork
!        logical(lk) :: overwrite_a_
!        logical, pointer :: bwork(:)
!        integer(ilp), allocatable :: iwork(:)
!        complex(dp), pointer :: amat(:,:), tau(:), work(:)
!        complex(dp) :: rwork_dummy(1)  ! Dummy for real/complex cases
!        complex(dp), allocatable :: tmat(:,:), zmat(:,:)
!        character :: jobz

!        ! Problem size
!        m = size(a, 1, kind=ilp)
!        lda = size(a, 1, kind=ilp)
!
!        ! Validate dimensions
!        if (size(a, 1, kind=ilp) /= size(a, 2, kind=ilp)) then
!            err0 = linalg_state(this, LINALG_VALUE_ERROR, 'Matrix A must be square: a=', [size(a,1), size(a,2)])
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Set default values
!        overwrite_a_ = .false._lk
!        if (present(overwrite_a)) overwrite_a_ = overwrite_a
!
!        ! Job type based on sorting request
!        if (present(sort) .and. sort) then
!            jobz = 'S'  ! Compute and sort eigenvalues
!        else
!            jobz = 'N'  ! Compute Schur decomposition without sorting
!        end if
!
!        ! Prepare storage
!        allocate(tmat(m, m), source=0.0_dp)
!        allocate(zmat(m, m), source=0.0_dp)
!
!        if (overwrite_a_) then
!            amat => a
!        else
!            allocate(amat(m, m), source=a)
!        end if
!
!        ! Allocate workspace
!        liwork = -1
!        if (present(lwork)) then
!            allocate(work(lwork))
!        else
!            allocate(work(1))  ! Temporary workspace for querying size
!        end if
!
!        ! Workspace query
!        call  zgees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, .false._lk, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Update workspace size
!        if (.not.present(lwork)) then
!            liwork = ceiling(real(work(1), kind=dp), kind=ilp)
!            deallocate(work)
!            allocate(work(liwork))
!        end if
!
!        ! Compute Schur decomposition
!        call  zgees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, present(sort) .and. sort, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Output results
!        t = amat
!        z = zmat

    end subroutine stdlib_linalg_z_schur

    subroutine get_schur_w_workspace(a,lwork,err)
        !> Input matrix a[m,m]
        complex(qp),intent(in),target :: a(:,:)
        !> Minimum workspace size for the decomposition operation
        integer(ilp),intent(out) :: lwork
        !> State return flag. Returns an error if the query failed
        type(linalg_state),optional,intent(out) :: err
        
        integer(ilp) :: m,n,sdim,info
        character :: jobvs,sort
        logical(lk) :: bwork_dummy(1)
        complex(qp),pointer :: amat(:,:)
        real(qp) :: rwork_dummy(1)
        complex(qp) :: wr_dummy(1),wi_dummy(1),vs_dummy(1,1),work_dummy(1)
        type(linalg_state) :: err0
        
        !> Initialize problem
        lwork = -1_ilp
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        
        !> Create a dummy intent(inout) argument
        amat => a
        
        !> Select dummy task
        jobvs = gees_vectors(.true.)
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        !> Get Schur workspace
        call gees(jobvs,sort,do_not_select,n,amat,m,sdim,wr_dummy, &
                  vs_dummy,m,work_dummy,lwork,rwork_dummy,bwork_dummy,info)
        if (info == 0) lwork = int(work_dummy(1),kind=ilp)
        call handle_gees_info(info,m,n,m,err0)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Interface to dummy select routine
            pure logical(lk) function do_not_select(alpha)
                complex(qp),intent(in) :: alpha
                do_not_select = .false.
            end function do_not_select
        
    end subroutine get_schur_w_workspace
    
    ! Schur decomposition subroutine
    pure subroutine stdlib_linalg_w_schur(a,t,z,storage,err)
        !> Input matrix a[m,m]
        complex(qp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(qp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(qp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(qp),intent(out),optional,target :: storage(:)
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
!        type(linalg_state) :: err0
!        integer(ilp) :: m, lda, info, liwork
!        logical(lk) :: overwrite_a_
!        logical, pointer :: bwork(:)
!        integer(ilp), allocatable :: iwork(:)
!        complex(qp), pointer :: amat(:,:), tau(:), work(:)
!        complex(qp) :: rwork_dummy(1)  ! Dummy for real/complex cases
!        complex(qp), allocatable :: tmat(:,:), zmat(:,:)
!        character :: jobz

!        ! Problem size
!        m = size(a, 1, kind=ilp)
!        lda = size(a, 1, kind=ilp)
!
!        ! Validate dimensions
!        if (size(a, 1, kind=ilp) /= size(a, 2, kind=ilp)) then
!            err0 = linalg_state(this, LINALG_VALUE_ERROR, 'Matrix A must be square: a=', [size(a,1), size(a,2)])
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Set default values
!        overwrite_a_ = .false._lk
!        if (present(overwrite_a)) overwrite_a_ = overwrite_a
!
!        ! Job type based on sorting request
!        if (present(sort) .and. sort) then
!            jobz = 'S'  ! Compute and sort eigenvalues
!        else
!            jobz = 'N'  ! Compute Schur decomposition without sorting
!        end if
!
!        ! Prepare storage
!        allocate(tmat(m, m), source=0.0_qp)
!        allocate(zmat(m, m), source=0.0_qp)
!
!        if (overwrite_a_) then
!            amat => a
!        else
!            allocate(amat(m, m), source=a)
!        end if
!
!        ! Allocate workspace
!        liwork = -1
!        if (present(lwork)) then
!            allocate(work(lwork))
!        else
!            allocate(work(1))  ! Temporary workspace for querying size
!        end if
!
!        ! Workspace query
!        call  zgees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, .false._lk, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Update workspace size
!        if (.not.present(lwork)) then
!            liwork = ceiling(real(work(1), kind=qp), kind=ilp)
!            deallocate(work)
!            allocate(work(liwork))
!        end if
!
!        ! Compute Schur decomposition
!        call  zgees  &
!            (jobz, 'N', nullify(bwork), m, amat, lda, tau, zmat, lda, work, liwork, iwork, info)
!        call handle_gees_info(info, m, present(sort) .and. sort, err0)
!        if (err0%error()) then
!            call linalg_error_handling(err0, err)
!            return
!        end if
!
!        ! Output results
!        t = amat
!        z = zmat

    end subroutine stdlib_linalg_w_schur

end module stdlib_linalg_schur

