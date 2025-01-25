module la_schur
    use la_constants
    use la_lapack,only:gees
    use la_state,only:linalg_state,linalg_error_handling,LINALG_ERROR, &
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
    !! ([Specification](../page/specs/la_linalg.html#schur-compute-the-schur-decomposition-of-a-matrix))
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
      module procedure la_s_schur
      module procedure la_real_eig_s_schur
      module procedure la_d_schur
      module procedure la_real_eig_d_schur
      module procedure la_q_schur
      module procedure la_real_eig_q_schur
      module procedure la_c_schur
      module procedure la_real_eig_c_schur
      module procedure la_z_schur
      module procedure la_real_eig_z_schur
      module procedure la_w_schur
      module procedure la_real_eig_w_schur
  end interface schur

  ! Return the working array space required by the Schur decomposition solver
  interface schur_space
    !! version: experimental
    !!
    !! Computes the working array space required by the Schur decomposition solver
    !! ([Specification](../page/specs/la_linalg.html#schur-space-compute-internal-working-space-requirements-for-the-schur-decomposition))
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
            elseif (info == n) then
                err = linalg_state(this,LINALG_ERROR,'Convergence failure: no converged eigenvalues')
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
        if (info == 0) lwork = nint(real(work_dummy(1),kind=sp),kind=ilp)
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
    subroutine la_s_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        real(sp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(sp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(sp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Output eigenvalues that appear on the diagonal of T
        complex(sp),optional,intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(sp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
        integer(ilp) :: m,n,mt,nt,ldvs,nvs,lde,lwork,sdim,info
        logical(lk) :: overwrite_a_
        logical(lk),target :: bwork_dummy(1),local_eigs
        logical(lk),pointer :: bwork(:)
        real(sp),allocatable :: rwork(:)
        real(sp),target :: vs_dummy(1,1)
        real(sp),pointer :: vs(:,:),work(:),eigs(:),eigi(:)
        character :: jobvs,sort
        type(linalg_state_type) :: err0
        
        ! Problem size
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        mt = size(t,1,kind=ilp)
        nt = size(t,2,kind=ilp)
        
        ! Validate dimensions
        if (m /= n .or. m <= 0 .or. n <= 0) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix A must be square: size(a)=', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        if (mt /= nt .or. mt /= n .or. nt /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix T must be square: size(T)=', [mt,nt], &
                                                          'should be', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        
        !> Copy data into the output array
        t = a
        
        ! Can A be overwritten? By default, do not overwrite
        overwrite_a_ = .false._lk
        if (present(overwrite_a)) overwrite_a_ = overwrite_a .and. n >= 2
        
        !> Schur vectors
        jobvs = gees_vectors(present(z))
        if (present(z)) then
            vs => z
            
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
            
            if (ldvs < n .or. nvs /= n) then
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Schur vectors size=', [ldvs,nvs], &
                                                              'should be n=',n)
                call linalg_error_handling(err0,err)
                return
            end if
            
        else
            vs => vs_dummy
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
        end if
        
        !> User or self-allocated storage
        if (present(storage)) then
            
            work => storage
            lwork = size(work,1,kind=ilp)
            
        else
            
            ! Query optimal workspace
            call get_schur_s_workspace(a,lwork,err0)
            
            if (err0%error()) then
                call linalg_error_handling(err0,err)
                return
            else
                allocate (work(lwork))
            end if
            
        end if
        
        !> SORTING: no sorting options are currently supported
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        if (sort /= GEES_NOT) then
            
           allocate (bwork(n),source=.false.)
        
        else
            
           bwork => bwork_dummy
            
        end if
        
        !> User or self-allocated eigenvalue storage
        if (present(eigvals)) then
            lde = size(eigvals,1,kind=ilp)
            local_eigs = .true.
        else
            local_eigs = .true.
            lde = n
        end if

        if (local_eigs) then
            ! Use A storage if possible
            if (overwrite_a_) then
                eigs => a(:,1)
                eigi => a(:,2)
            else
                allocate (eigs(n),eigi(n))
            end if
        end if
        
        if (lde < n) then
            
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, &
                                           'Insufficient eigenvalue array size=',lde, &
                                           'should be >=',n)
        
        else
            
            ! Compute Schur decomposition
            call gees(jobvs,sort,eig_select,nt,t,mt,sdim,eigs,eigi, &
                      vs,ldvs,work,lwork,bwork,info)
            call handle_gees_info(info,m,n,m,err0)
            
        end if

        eigenvalue_output: if (local_eigs) then
           ! Build complex eigenvalues
           if (present(eigvals)) eigvals = cmplx(eigs,eigi,kind=sp)
           if (.not. overwrite_a_) deallocate (eigs,eigi)
        end if eigenvalue_output
        if (.not. present(storage)) deallocate (work)
        if (sort /= GEES_NOT) deallocate (bwork)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Dummy select routine: currently, no sorting options are offered
            pure logical(lk) function eig_select(alphar,alphai)
                real(sp),intent(in) :: alphar,alphai
                complex(sp) :: alpha
                alpha = cmplx(alphar,alphai,kind=sp)
                eig_select = .false.
            end function eig_select

    end subroutine la_s_schur

    ! Schur decomposition subroutine: real eigenvalue interface
    module subroutine la_real_eig_s_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        real(sp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(sp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(sp),optional,intent(out),contiguous,target :: z(:,:)
        !> Output eigenvalues that appear on the diagonal of T
        real(sp),intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(sp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err
        
        type(linalg_state) :: err0
        integer(ilp) :: n
        complex(sp),allocatable :: ceigvals(:)
        real(sp),parameter :: rtol = epsilon(0.0_sp)
        real(sp),parameter :: atol = tiny(0.0_sp)
          
        n = size(eigvals,dim=1,kind=ilp)
        allocate (ceigvals(n))
          
        !> Compute Schur decomposition with complex eigenvalues
        call la_s_schur(a,t,z,ceigvals,overwrite_a,storage,err0)
          
        ! Check that no eigenvalues have meaningful imaginary part
        if (err0%ok() .and. any(aimag(ceigvals) > atol + rtol*abs(abs(ceigvals)))) then
           err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                              'complex eigenvalues detected: max(imag(lambda))=',maxval(aimag(ceigvals)))
        end if
          
        ! Return real components only
        eigvals(:n) = real(ceigvals,kind=sp)
          
        call linalg_error_handling(err0,err)
        
    end subroutine la_real_eig_s_schur

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
        if (info == 0) lwork = nint(real(work_dummy(1),kind=dp),kind=ilp)
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
    subroutine la_d_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        real(dp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(dp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(dp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Output eigenvalues that appear on the diagonal of T
        complex(dp),optional,intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(dp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
        integer(ilp) :: m,n,mt,nt,ldvs,nvs,lde,lwork,sdim,info
        logical(lk) :: overwrite_a_
        logical(lk),target :: bwork_dummy(1),local_eigs
        logical(lk),pointer :: bwork(:)
        real(dp),allocatable :: rwork(:)
        real(dp),target :: vs_dummy(1,1)
        real(dp),pointer :: vs(:,:),work(:),eigs(:),eigi(:)
        character :: jobvs,sort
        type(linalg_state_type) :: err0
        
        ! Problem size
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        mt = size(t,1,kind=ilp)
        nt = size(t,2,kind=ilp)
        
        ! Validate dimensions
        if (m /= n .or. m <= 0 .or. n <= 0) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix A must be square: size(a)=', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        if (mt /= nt .or. mt /= n .or. nt /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix T must be square: size(T)=', [mt,nt], &
                                                          'should be', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        
        !> Copy data into the output array
        t = a
        
        ! Can A be overwritten? By default, do not overwrite
        overwrite_a_ = .false._lk
        if (present(overwrite_a)) overwrite_a_ = overwrite_a .and. n >= 2
        
        !> Schur vectors
        jobvs = gees_vectors(present(z))
        if (present(z)) then
            vs => z
            
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
            
            if (ldvs < n .or. nvs /= n) then
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Schur vectors size=', [ldvs,nvs], &
                                                              'should be n=',n)
                call linalg_error_handling(err0,err)
                return
            end if
            
        else
            vs => vs_dummy
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
        end if
        
        !> User or self-allocated storage
        if (present(storage)) then
            
            work => storage
            lwork = size(work,1,kind=ilp)
            
        else
            
            ! Query optimal workspace
            call get_schur_d_workspace(a,lwork,err0)
            
            if (err0%error()) then
                call linalg_error_handling(err0,err)
                return
            else
                allocate (work(lwork))
            end if
            
        end if
        
        !> SORTING: no sorting options are currently supported
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        if (sort /= GEES_NOT) then
            
           allocate (bwork(n),source=.false.)
        
        else
            
           bwork => bwork_dummy
            
        end if
        
        !> User or self-allocated eigenvalue storage
        if (present(eigvals)) then
            lde = size(eigvals,1,kind=ilp)
            local_eigs = .true.
        else
            local_eigs = .true.
            lde = n
        end if

        if (local_eigs) then
            ! Use A storage if possible
            if (overwrite_a_) then
                eigs => a(:,1)
                eigi => a(:,2)
            else
                allocate (eigs(n),eigi(n))
            end if
        end if
        
        if (lde < n) then
            
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, &
                                           'Insufficient eigenvalue array size=',lde, &
                                           'should be >=',n)
        
        else
            
            ! Compute Schur decomposition
            call gees(jobvs,sort,eig_select,nt,t,mt,sdim,eigs,eigi, &
                      vs,ldvs,work,lwork,bwork,info)
            call handle_gees_info(info,m,n,m,err0)
            
        end if

        eigenvalue_output: if (local_eigs) then
           ! Build complex eigenvalues
           if (present(eigvals)) eigvals = cmplx(eigs,eigi,kind=dp)
           if (.not. overwrite_a_) deallocate (eigs,eigi)
        end if eigenvalue_output
        if (.not. present(storage)) deallocate (work)
        if (sort /= GEES_NOT) deallocate (bwork)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Dummy select routine: currently, no sorting options are offered
            pure logical(lk) function eig_select(alphar,alphai)
                real(dp),intent(in) :: alphar,alphai
                complex(dp) :: alpha
                alpha = cmplx(alphar,alphai,kind=dp)
                eig_select = .false.
            end function eig_select

    end subroutine la_d_schur

    ! Schur decomposition subroutine: real eigenvalue interface
    module subroutine la_real_eig_d_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        real(dp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(dp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(dp),optional,intent(out),contiguous,target :: z(:,:)
        !> Output eigenvalues that appear on the diagonal of T
        real(dp),intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(dp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err
        
        type(linalg_state) :: err0
        integer(ilp) :: n
        complex(dp),allocatable :: ceigvals(:)
        real(dp),parameter :: rtol = epsilon(0.0_dp)
        real(dp),parameter :: atol = tiny(0.0_dp)
          
        n = size(eigvals,dim=1,kind=ilp)
        allocate (ceigvals(n))
          
        !> Compute Schur decomposition with complex eigenvalues
        call la_d_schur(a,t,z,ceigvals,overwrite_a,storage,err0)
          
        ! Check that no eigenvalues have meaningful imaginary part
        if (err0%ok() .and. any(aimag(ceigvals) > atol + rtol*abs(abs(ceigvals)))) then
           err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                              'complex eigenvalues detected: max(imag(lambda))=',maxval(aimag(ceigvals)))
        end if
          
        ! Return real components only
        eigvals(:n) = real(ceigvals,kind=dp)
          
        call linalg_error_handling(err0,err)
        
    end subroutine la_real_eig_d_schur

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
        if (info == 0) lwork = nint(real(work_dummy(1),kind=qp),kind=ilp)
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
    subroutine la_q_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        real(qp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(qp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(qp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Output eigenvalues that appear on the diagonal of T
        complex(qp),optional,intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(qp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
        integer(ilp) :: m,n,mt,nt,ldvs,nvs,lde,lwork,sdim,info
        logical(lk) :: overwrite_a_
        logical(lk),target :: bwork_dummy(1),local_eigs
        logical(lk),pointer :: bwork(:)
        real(qp),allocatable :: rwork(:)
        real(qp),target :: vs_dummy(1,1)
        real(qp),pointer :: vs(:,:),work(:),eigs(:),eigi(:)
        character :: jobvs,sort
        type(linalg_state_type) :: err0
        
        ! Problem size
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        mt = size(t,1,kind=ilp)
        nt = size(t,2,kind=ilp)
        
        ! Validate dimensions
        if (m /= n .or. m <= 0 .or. n <= 0) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix A must be square: size(a)=', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        if (mt /= nt .or. mt /= n .or. nt /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix T must be square: size(T)=', [mt,nt], &
                                                          'should be', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        
        !> Copy data into the output array
        t = a
        
        ! Can A be overwritten? By default, do not overwrite
        overwrite_a_ = .false._lk
        if (present(overwrite_a)) overwrite_a_ = overwrite_a .and. n >= 2
        
        !> Schur vectors
        jobvs = gees_vectors(present(z))
        if (present(z)) then
            vs => z
            
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
            
            if (ldvs < n .or. nvs /= n) then
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Schur vectors size=', [ldvs,nvs], &
                                                              'should be n=',n)
                call linalg_error_handling(err0,err)
                return
            end if
            
        else
            vs => vs_dummy
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
        end if
        
        !> User or self-allocated storage
        if (present(storage)) then
            
            work => storage
            lwork = size(work,1,kind=ilp)
            
        else
            
            ! Query optimal workspace
            call get_schur_q_workspace(a,lwork,err0)
            
            if (err0%error()) then
                call linalg_error_handling(err0,err)
                return
            else
                allocate (work(lwork))
            end if
            
        end if
        
        !> SORTING: no sorting options are currently supported
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        if (sort /= GEES_NOT) then
            
           allocate (bwork(n),source=.false.)
        
        else
            
           bwork => bwork_dummy
            
        end if
        
        !> User or self-allocated eigenvalue storage
        if (present(eigvals)) then
            lde = size(eigvals,1,kind=ilp)
            local_eigs = .true.
        else
            local_eigs = .true.
            lde = n
        end if

        if (local_eigs) then
            ! Use A storage if possible
            if (overwrite_a_) then
                eigs => a(:,1)
                eigi => a(:,2)
            else
                allocate (eigs(n),eigi(n))
            end if
        end if
        
        if (lde < n) then
            
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, &
                                           'Insufficient eigenvalue array size=',lde, &
                                           'should be >=',n)
        
        else
            
            ! Compute Schur decomposition
            call gees(jobvs,sort,eig_select,nt,t,mt,sdim,eigs,eigi, &
                      vs,ldvs,work,lwork,bwork,info)
            call handle_gees_info(info,m,n,m,err0)
            
        end if

        eigenvalue_output: if (local_eigs) then
           ! Build complex eigenvalues
           if (present(eigvals)) eigvals = cmplx(eigs,eigi,kind=qp)
           if (.not. overwrite_a_) deallocate (eigs,eigi)
        end if eigenvalue_output
        if (.not. present(storage)) deallocate (work)
        if (sort /= GEES_NOT) deallocate (bwork)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Dummy select routine: currently, no sorting options are offered
            pure logical(lk) function eig_select(alphar,alphai)
                real(qp),intent(in) :: alphar,alphai
                complex(qp) :: alpha
                alpha = cmplx(alphar,alphai,kind=qp)
                eig_select = .false.
            end function eig_select

    end subroutine la_q_schur

    ! Schur decomposition subroutine: real eigenvalue interface
    module subroutine la_real_eig_q_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        real(qp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        real(qp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        real(qp),optional,intent(out),contiguous,target :: z(:,:)
        !> Output eigenvalues that appear on the diagonal of T
        real(qp),intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        real(qp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err
        
        type(linalg_state) :: err0
        integer(ilp) :: n
        complex(qp),allocatable :: ceigvals(:)
        real(qp),parameter :: rtol = epsilon(0.0_qp)
        real(qp),parameter :: atol = tiny(0.0_qp)
          
        n = size(eigvals,dim=1,kind=ilp)
        allocate (ceigvals(n))
          
        !> Compute Schur decomposition with complex eigenvalues
        call la_q_schur(a,t,z,ceigvals,overwrite_a,storage,err0)
          
        ! Check that no eigenvalues have meaningful imaginary part
        if (err0%ok() .and. any(aimag(ceigvals) > atol + rtol*abs(abs(ceigvals)))) then
           err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                              'complex eigenvalues detected: max(imag(lambda))=',maxval(aimag(ceigvals)))
        end if
          
        ! Return real components only
        eigvals(:n) = real(ceigvals,kind=qp)
          
        call linalg_error_handling(err0,err)
        
    end subroutine la_real_eig_q_schur

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
        if (info == 0) lwork = nint(real(work_dummy(1),kind=sp),kind=ilp)
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
    subroutine la_c_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        complex(sp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(sp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(sp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Output eigenvalues that appear on the diagonal of T
        complex(sp),optional,intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(sp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
        integer(ilp) :: m,n,mt,nt,ldvs,nvs,lde,lwork,sdim,info
        logical(lk) :: overwrite_a_
        logical(lk),target :: bwork_dummy(1),local_eigs
        logical(lk),pointer :: bwork(:)
        real(sp),allocatable :: rwork(:)
        complex(sp),target :: vs_dummy(1,1)
        complex(sp),pointer :: vs(:,:),work(:),eigs(:)
        character :: jobvs,sort
        type(linalg_state_type) :: err0
        
        ! Problem size
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        mt = size(t,1,kind=ilp)
        nt = size(t,2,kind=ilp)
        
        ! Validate dimensions
        if (m /= n .or. m <= 0 .or. n <= 0) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix A must be square: size(a)=', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        if (mt /= nt .or. mt /= n .or. nt /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix T must be square: size(T)=', [mt,nt], &
                                                          'should be', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        
        !> Copy data into the output array
        t = a
        
        ! Can A be overwritten? By default, do not overwrite
        overwrite_a_ = .false._lk
        if (present(overwrite_a)) overwrite_a_ = overwrite_a .and. n >= 2
        
        !> Schur vectors
        jobvs = gees_vectors(present(z))
        if (present(z)) then
            vs => z
            
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
            
            if (ldvs < n .or. nvs /= n) then
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Schur vectors size=', [ldvs,nvs], &
                                                              'should be n=',n)
                call linalg_error_handling(err0,err)
                return
            end if
            
        else
            vs => vs_dummy
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
        end if
        
        !> User or self-allocated storage
        if (present(storage)) then
            
            work => storage
            lwork = size(work,1,kind=ilp)
            
        else
            
            ! Query optimal workspace
            call get_schur_c_workspace(a,lwork,err0)
            
            if (err0%error()) then
                call linalg_error_handling(err0,err)
                return
            else
                allocate (work(lwork))
            end if
            
        end if
        
        !> SORTING: no sorting options are currently supported
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        if (sort /= GEES_NOT) then
            
           allocate (bwork(n),source=.false.)
        
        else
            
           bwork => bwork_dummy
            
        end if
        
        !> User or self-allocated eigenvalue storage
        if (present(eigvals)) then
            lde = size(eigvals,1,kind=ilp)
            eigs => eigvals
            local_eigs = .false.
        else
            local_eigs = .true.
            lde = n
        end if

        if (local_eigs) then
            ! Use A storage if possible
            if (overwrite_a_) then
                eigs => a(:,1)
            else
                allocate (eigs(n))
            end if
        end if
        
        allocate (rwork(n))

        if (lde < n) then
            
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, &
                                           'Insufficient eigenvalue array size=',lde, &
                                           'should be >=',n)
        
        else
            
            ! Compute Schur decomposition
            call gees(jobvs,sort,eig_select,nt,t,mt,sdim,eigs, &
                      vs,ldvs,work,lwork,rwork,bwork,info)
            call handle_gees_info(info,m,n,m,err0)
            
        end if

        eigenvalue_output: if (local_eigs) then
           if (.not. overwrite_a_) deallocate (eigs)
        end if eigenvalue_output
        if (.not. present(storage)) deallocate (work)
        if (sort /= GEES_NOT) deallocate (bwork)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Dummy select routine: currently, no sorting options are offered
            pure logical(lk) function eig_select(alpha)
                complex(sp),intent(in) :: alpha
                eig_select = .false.
            end function eig_select

    end subroutine la_c_schur

    ! Schur decomposition subroutine: real eigenvalue interface
    module subroutine la_real_eig_c_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        complex(sp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(sp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(sp),optional,intent(out),contiguous,target :: z(:,:)
        !> Output eigenvalues that appear on the diagonal of T
        real(sp),intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(sp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err
        
        type(linalg_state) :: err0
        integer(ilp) :: n
        complex(sp),allocatable :: ceigvals(:)
        real(sp),parameter :: rtol = epsilon(0.0_sp)
        real(sp),parameter :: atol = tiny(0.0_sp)
          
        n = size(eigvals,dim=1,kind=ilp)
        allocate (ceigvals(n))
          
        !> Compute Schur decomposition with complex eigenvalues
        call la_c_schur(a,t,z,ceigvals,overwrite_a,storage,err0)
          
        ! Check that no eigenvalues have meaningful imaginary part
        if (err0%ok() .and. any(aimag(ceigvals) > atol + rtol*abs(abs(ceigvals)))) then
           err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                              'complex eigenvalues detected: max(imag(lambda))=',maxval(aimag(ceigvals)))
        end if
          
        ! Return real components only
        eigvals(:n) = real(ceigvals,kind=sp)
          
        call linalg_error_handling(err0,err)
        
    end subroutine la_real_eig_c_schur

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
        if (info == 0) lwork = nint(real(work_dummy(1),kind=dp),kind=ilp)
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
    subroutine la_z_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        complex(dp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(dp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(dp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Output eigenvalues that appear on the diagonal of T
        complex(dp),optional,intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(dp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
        integer(ilp) :: m,n,mt,nt,ldvs,nvs,lde,lwork,sdim,info
        logical(lk) :: overwrite_a_
        logical(lk),target :: bwork_dummy(1),local_eigs
        logical(lk),pointer :: bwork(:)
        real(dp),allocatable :: rwork(:)
        complex(dp),target :: vs_dummy(1,1)
        complex(dp),pointer :: vs(:,:),work(:),eigs(:)
        character :: jobvs,sort
        type(linalg_state_type) :: err0
        
        ! Problem size
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        mt = size(t,1,kind=ilp)
        nt = size(t,2,kind=ilp)
        
        ! Validate dimensions
        if (m /= n .or. m <= 0 .or. n <= 0) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix A must be square: size(a)=', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        if (mt /= nt .or. mt /= n .or. nt /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix T must be square: size(T)=', [mt,nt], &
                                                          'should be', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        
        !> Copy data into the output array
        t = a
        
        ! Can A be overwritten? By default, do not overwrite
        overwrite_a_ = .false._lk
        if (present(overwrite_a)) overwrite_a_ = overwrite_a .and. n >= 2
        
        !> Schur vectors
        jobvs = gees_vectors(present(z))
        if (present(z)) then
            vs => z
            
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
            
            if (ldvs < n .or. nvs /= n) then
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Schur vectors size=', [ldvs,nvs], &
                                                              'should be n=',n)
                call linalg_error_handling(err0,err)
                return
            end if
            
        else
            vs => vs_dummy
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
        end if
        
        !> User or self-allocated storage
        if (present(storage)) then
            
            work => storage
            lwork = size(work,1,kind=ilp)
            
        else
            
            ! Query optimal workspace
            call get_schur_z_workspace(a,lwork,err0)
            
            if (err0%error()) then
                call linalg_error_handling(err0,err)
                return
            else
                allocate (work(lwork))
            end if
            
        end if
        
        !> SORTING: no sorting options are currently supported
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        if (sort /= GEES_NOT) then
            
           allocate (bwork(n),source=.false.)
        
        else
            
           bwork => bwork_dummy
            
        end if
        
        !> User or self-allocated eigenvalue storage
        if (present(eigvals)) then
            lde = size(eigvals,1,kind=ilp)
            eigs => eigvals
            local_eigs = .false.
        else
            local_eigs = .true.
            lde = n
        end if

        if (local_eigs) then
            ! Use A storage if possible
            if (overwrite_a_) then
                eigs => a(:,1)
            else
                allocate (eigs(n))
            end if
        end if
        
        allocate (rwork(n))

        if (lde < n) then
            
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, &
                                           'Insufficient eigenvalue array size=',lde, &
                                           'should be >=',n)
        
        else
            
            ! Compute Schur decomposition
            call gees(jobvs,sort,eig_select,nt,t,mt,sdim,eigs, &
                      vs,ldvs,work,lwork,rwork,bwork,info)
            call handle_gees_info(info,m,n,m,err0)
            
        end if

        eigenvalue_output: if (local_eigs) then
           if (.not. overwrite_a_) deallocate (eigs)
        end if eigenvalue_output
        if (.not. present(storage)) deallocate (work)
        if (sort /= GEES_NOT) deallocate (bwork)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Dummy select routine: currently, no sorting options are offered
            pure logical(lk) function eig_select(alpha)
                complex(dp),intent(in) :: alpha
                eig_select = .false.
            end function eig_select

    end subroutine la_z_schur

    ! Schur decomposition subroutine: real eigenvalue interface
    module subroutine la_real_eig_z_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        complex(dp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(dp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(dp),optional,intent(out),contiguous,target :: z(:,:)
        !> Output eigenvalues that appear on the diagonal of T
        real(dp),intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(dp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err
        
        type(linalg_state) :: err0
        integer(ilp) :: n
        complex(dp),allocatable :: ceigvals(:)
        real(dp),parameter :: rtol = epsilon(0.0_dp)
        real(dp),parameter :: atol = tiny(0.0_dp)
          
        n = size(eigvals,dim=1,kind=ilp)
        allocate (ceigvals(n))
          
        !> Compute Schur decomposition with complex eigenvalues
        call la_z_schur(a,t,z,ceigvals,overwrite_a,storage,err0)
          
        ! Check that no eigenvalues have meaningful imaginary part
        if (err0%ok() .and. any(aimag(ceigvals) > atol + rtol*abs(abs(ceigvals)))) then
           err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                              'complex eigenvalues detected: max(imag(lambda))=',maxval(aimag(ceigvals)))
        end if
          
        ! Return real components only
        eigvals(:n) = real(ceigvals,kind=dp)
          
        call linalg_error_handling(err0,err)
        
    end subroutine la_real_eig_z_schur

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
        if (info == 0) lwork = nint(real(work_dummy(1),kind=qp),kind=ilp)
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
    subroutine la_w_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        complex(qp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(qp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(qp),optional,intent(out),contiguous,target :: z(:,:)
        !> [optional] Output eigenvalues that appear on the diagonal of T
        complex(qp),optional,intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(qp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err

        ! Local variables
        integer(ilp) :: m,n,mt,nt,ldvs,nvs,lde,lwork,sdim,info
        logical(lk) :: overwrite_a_
        logical(lk),target :: bwork_dummy(1),local_eigs
        logical(lk),pointer :: bwork(:)
        real(qp),allocatable :: rwork(:)
        complex(qp),target :: vs_dummy(1,1)
        complex(qp),pointer :: vs(:,:),work(:),eigs(:)
        character :: jobvs,sort
        type(linalg_state_type) :: err0
        
        ! Problem size
        m = size(a,1,kind=ilp)
        n = size(a,2,kind=ilp)
        mt = size(t,1,kind=ilp)
        nt = size(t,2,kind=ilp)
        
        ! Validate dimensions
        if (m /= n .or. m <= 0 .or. n <= 0) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix A must be square: size(a)=', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        if (mt /= nt .or. mt /= n .or. nt /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Matrix T must be square: size(T)=', [mt,nt], &
                                                          'should be', [m,n])
            call linalg_error_handling(err0,err)
            return
        end if
        
        !> Copy data into the output array
        t = a
        
        ! Can A be overwritten? By default, do not overwrite
        overwrite_a_ = .false._lk
        if (present(overwrite_a)) overwrite_a_ = overwrite_a .and. n >= 2
        
        !> Schur vectors
        jobvs = gees_vectors(present(z))
        if (present(z)) then
            vs => z
            
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
            
            if (ldvs < n .or. nvs /= n) then
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Schur vectors size=', [ldvs,nvs], &
                                                              'should be n=',n)
                call linalg_error_handling(err0,err)
                return
            end if
            
        else
            vs => vs_dummy
            ldvs = size(vs,1,kind=ilp)
            nvs = size(vs,2,kind=ilp)
        end if
        
        !> User or self-allocated storage
        if (present(storage)) then
            
            work => storage
            lwork = size(work,1,kind=ilp)
            
        else
            
            ! Query optimal workspace
            call get_schur_w_workspace(a,lwork,err0)
            
            if (err0%error()) then
                call linalg_error_handling(err0,err)
                return
            else
                allocate (work(lwork))
            end if
            
        end if
        
        !> SORTING: no sorting options are currently supported
        sort = gees_sort_eigs(.false.)
        sdim = 0_ilp
        
        if (sort /= GEES_NOT) then
            
           allocate (bwork(n),source=.false.)
        
        else
            
           bwork => bwork_dummy
            
        end if
        
        !> User or self-allocated eigenvalue storage
        if (present(eigvals)) then
            lde = size(eigvals,1,kind=ilp)
            eigs => eigvals
            local_eigs = .false.
        else
            local_eigs = .true.
            lde = n
        end if

        if (local_eigs) then
            ! Use A storage if possible
            if (overwrite_a_) then
                eigs => a(:,1)
            else
                allocate (eigs(n))
            end if
        end if
        
        allocate (rwork(n))

        if (lde < n) then
            
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, &
                                           'Insufficient eigenvalue array size=',lde, &
                                           'should be >=',n)
        
        else
            
            ! Compute Schur decomposition
            call gees(jobvs,sort,eig_select,nt,t,mt,sdim,eigs, &
                      vs,ldvs,work,lwork,rwork,bwork,info)
            call handle_gees_info(info,m,n,m,err0)
            
        end if

        eigenvalue_output: if (local_eigs) then
           if (.not. overwrite_a_) deallocate (eigs)
        end if eigenvalue_output
        if (.not. present(storage)) deallocate (work)
        if (sort /= GEES_NOT) deallocate (bwork)
        call linalg_error_handling(err0,err)
        
        contains
        
            ! Dummy select routine: currently, no sorting options are offered
            pure logical(lk) function eig_select(alpha)
                complex(qp),intent(in) :: alpha
                eig_select = .false.
            end function eig_select

    end subroutine la_w_schur

    ! Schur decomposition subroutine: real eigenvalue interface
    module subroutine la_real_eig_w_schur(a,t,z,eigvals,overwrite_a,storage,err)
        !> Input matrix a[m,m]
        complex(qp),intent(inout),target :: a(:,:)
        !> Schur form of A: upper-triangular or quasi-upper-triangular matrix T
        complex(qp),intent(out),contiguous,target :: t(:,:)
        !> Unitary/orthonormal transformation matrix Z
        complex(qp),optional,intent(out),contiguous,target :: z(:,:)
        !> Output eigenvalues that appear on the diagonal of T
        real(qp),intent(out),contiguous,target :: eigvals(:)
        !> [optional] Provide pre-allocated workspace, size to be checked with schur_space
        complex(qp),optional,intent(inout),target :: storage(:)
        !> [optional] Can A data be overwritten and destroyed?
        logical(lk),optional,intent(in) :: overwrite_a
        !> [optional] State return flag. On error if not requested, the code will stop
        type(linalg_state),optional,intent(out) :: err
        
        type(linalg_state) :: err0
        integer(ilp) :: n
        complex(qp),allocatable :: ceigvals(:)
        real(qp),parameter :: rtol = epsilon(0.0_qp)
        real(qp),parameter :: atol = tiny(0.0_qp)
          
        n = size(eigvals,dim=1,kind=ilp)
        allocate (ceigvals(n))
          
        !> Compute Schur decomposition with complex eigenvalues
        call la_w_schur(a,t,z,ceigvals,overwrite_a,storage,err0)
          
        ! Check that no eigenvalues have meaningful imaginary part
        if (err0%ok() .and. any(aimag(ceigvals) > atol + rtol*abs(abs(ceigvals)))) then
           err0 = linalg_state(this,LINALG_VALUE_ERROR, &
                              'complex eigenvalues detected: max(imag(lambda))=',maxval(aimag(ceigvals)))
        end if
          
        ! Return real components only
        eigvals(:n) = real(ceigvals,kind=qp)
          
        call linalg_error_handling(err0,err)
        
    end subroutine la_real_eig_w_schur

end module la_schur

