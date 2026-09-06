!> Inverse of a square matrix
module la_inverse
     use la_constants
     use la_blas
     use la_lapack
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> @brief Compute the inverse of a square matrix.
     !!
     !! This function computes the inverse of a real or complex square matrix \f$ A \f$.
     !! The inverse is computed using an LU decomposition with partial pivoting.
     !!
     !! @param[in] A The input square matrix of size \f$ [n, n] \f$.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return The inverse matrix \f$ A^{-1} \f$ of size \f$ [n, n] \f$.
     !!
     !! @note This function relies on LAPACK's LU decomposition routines ([GETRF](@ref la_lapack::getrf)
     !!       and [GETRI](@ref la_lapack::getri)).
     !! @warning The matrix \f$ A \f$ must be non-singular. If it is singular or nearly singular,
     !!          the function will fail.
     !!
     public :: inv

     !> @brief Compute the inverse of a square matrix.
     !!
     !! This subroutine computes the inverse of a real or complex square matrix \f$ A \f$, either
     !! in-place or into a second matrix `inva` of the same shape.
     !! The inverse is computed using an LU decomposition with partial pivoting.
     !! Storage for the pivot indices may be provided, so that a repeated call performs no
     !! internal allocation of the pivot array.
     !!
     !! @param[in,out] A The input square matrix of size \f$ [n, n] \f$. In the in-place form it is
     !!                  replaced by its inverse \f$ A^{-1} \f$; in the split form it is read only.
     !! @param[out] inva (Split form only) The inverse matrix \f$ A^{-1} \f$, of size \f$ [n, n] \f$.
     !! @param[in,out] pivot (Optional) Storage array for the `n` diagonal pivot indices.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @note The in-place form is useful when memory efficiency is a priority, as it avoids additional allocations.
     !! @warning The matrix \f$ A \f$ must be non-singular. If it is singular or nearly singular,
     !!          the computation will fail.
     !!
     public :: invert

     !> @brief Compute the inverse of a square matrix using the `.inv.` operator.
     !!
     !! This operator computes the inverse of a real or complex square matrix \f$ A \f$ using
     !! an LU decomposition with partial pivoting.
     !!
     !! @param[in] A The input square matrix of size \f$ [n, n] \f$.
     !!
     !! @return The inverse matrix \f$ A^{-1} \f$ of size \f$ [n, n] \f$.
     !!
     !! @note This operator is a shorthand for calling `inv(A)`, allowing expressions such as:
     !!       \f$ X = .inv.A \f$
     !! @warning The matrix \f$ A \f$ must be non-singular. If it is singular or nearly singular,
     !!          the computation will fail.
     !!
     public :: operator(.inv.)

     ! Function interface
     interface inv
        module procedure la_inverse_s
        module procedure la_inverse_d
        module procedure la_inverse_q
        module procedure la_inverse_c
        module procedure la_inverse_z
        module procedure la_inverse_w
     end interface inv

     ! Subroutine interface: in-place and split factorization
     interface invert
        module procedure la_invert_s
        module procedure la_invert_split_s
        module procedure la_invert_d
        module procedure la_invert_split_d
        module procedure la_invert_q
        module procedure la_invert_split_q
        module procedure la_invert_c
        module procedure la_invert_split_c
        module procedure la_invert_z
        module procedure la_invert_split_z
        module procedure la_invert_w
        module procedure la_invert_split_w
     end interface invert

     ! Operator interface
     interface operator(.inv.)
        module procedure la_inverse_s_operator
        module procedure la_inverse_d_operator
        module procedure la_inverse_q_operator
        module procedure la_inverse_c_operator
        module procedure la_inverse_z_operator
        module procedure la_inverse_w_operator
     end interface operator(.inv.)

     contains

     ! Compute the in-place square matrix inverse of a
     subroutine la_invert_s(a,pivot,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout) :: a(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork,npiv
         integer(ilp),pointer :: ipiv(:)
         real(sp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n .or. npiv < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = la_ilaenv(1,'sgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
         if (.not. present(pivot)) deallocate (ipiv)
         call err0%handle(err)

     end subroutine la_invert_s

     ! Compute the square matrix inverse of a into a second matrix
     subroutine la_invert_split_s(a,inva,pivot,err)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Inverse matrix inva[n,n]
         real(sp),intent(out) :: inva(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: sa(2),sinva(2)
         character(*),parameter :: this = 'invert'

         sa = shape(a,kind=ilp)
         sinva = shape(inva,kind=ilp)

         if (any(sa /= sinva)) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',sa(1),',',sa(2),'],', &
                                                                       'inva=[',sinva(1),',',sinva(2),']')

         else

            !> Copy data in
            inva = a

            !> Compute matrix inverse
            call la_invert_s(inva,pivot=pivot,err=err0)

         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_invert_split_s

     ! Invert matrix in place
     function la_inverse_s(a,err) result(inva)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Output matrix inverse
         real(sp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call la_invert_s(inva,err=err)

     end function la_inverse_s

     ! Inverse matrix operator
     function la_inverse_s_operator(a) result(inva)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Result matrix
         real(sp),allocatable :: inva(:,:)

         type(la_state) :: err

         inva = la_inverse_s(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function la_inverse_s_operator

     ! Compute the in-place square matrix inverse of a
     subroutine la_invert_d(a,pivot,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout) :: a(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork,npiv
         integer(ilp),pointer :: ipiv(:)
         real(dp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n .or. npiv < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = la_ilaenv(1,'dgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
         if (.not. present(pivot)) deallocate (ipiv)
         call err0%handle(err)

     end subroutine la_invert_d

     ! Compute the square matrix inverse of a into a second matrix
     subroutine la_invert_split_d(a,inva,pivot,err)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Inverse matrix inva[n,n]
         real(dp),intent(out) :: inva(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: sa(2),sinva(2)
         character(*),parameter :: this = 'invert'

         sa = shape(a,kind=ilp)
         sinva = shape(inva,kind=ilp)

         if (any(sa /= sinva)) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',sa(1),',',sa(2),'],', &
                                                                       'inva=[',sinva(1),',',sinva(2),']')

         else

            !> Copy data in
            inva = a

            !> Compute matrix inverse
            call la_invert_d(inva,pivot=pivot,err=err0)

         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_invert_split_d

     ! Invert matrix in place
     function la_inverse_d(a,err) result(inva)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Output matrix inverse
         real(dp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call la_invert_d(inva,err=err)

     end function la_inverse_d

     ! Inverse matrix operator
     function la_inverse_d_operator(a) result(inva)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Result matrix
         real(dp),allocatable :: inva(:,:)

         type(la_state) :: err

         inva = la_inverse_d(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function la_inverse_d_operator

     ! Compute the in-place square matrix inverse of a
     subroutine la_invert_q(a,pivot,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout) :: a(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork,npiv
         integer(ilp),pointer :: ipiv(:)
         real(qp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n .or. npiv < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = la_ilaenv(1,'qgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
         if (.not. present(pivot)) deallocate (ipiv)
         call err0%handle(err)

     end subroutine la_invert_q

     ! Compute the square matrix inverse of a into a second matrix
     subroutine la_invert_split_q(a,inva,pivot,err)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Inverse matrix inva[n,n]
         real(qp),intent(out) :: inva(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: sa(2),sinva(2)
         character(*),parameter :: this = 'invert'

         sa = shape(a,kind=ilp)
         sinva = shape(inva,kind=ilp)

         if (any(sa /= sinva)) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',sa(1),',',sa(2),'],', &
                                                                       'inva=[',sinva(1),',',sinva(2),']')

         else

            !> Copy data in
            inva = a

            !> Compute matrix inverse
            call la_invert_q(inva,pivot=pivot,err=err0)

         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_invert_split_q

     ! Invert matrix in place
     function la_inverse_q(a,err) result(inva)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Output matrix inverse
         real(qp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call la_invert_q(inva,err=err)

     end function la_inverse_q

     ! Inverse matrix operator
     function la_inverse_q_operator(a) result(inva)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Result matrix
         real(qp),allocatable :: inva(:,:)

         type(la_state) :: err

         inva = la_inverse_q(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function la_inverse_q_operator

     ! Compute the in-place square matrix inverse of a
     subroutine la_invert_c(a,pivot,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout) :: a(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork,npiv
         integer(ilp),pointer :: ipiv(:)
         complex(sp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n .or. npiv < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = la_ilaenv(1,'cgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
         if (.not. present(pivot)) deallocate (ipiv)
         call err0%handle(err)

     end subroutine la_invert_c

     ! Compute the square matrix inverse of a into a second matrix
     subroutine la_invert_split_c(a,inva,pivot,err)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Inverse matrix inva[n,n]
         complex(sp),intent(out) :: inva(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: sa(2),sinva(2)
         character(*),parameter :: this = 'invert'

         sa = shape(a,kind=ilp)
         sinva = shape(inva,kind=ilp)

         if (any(sa /= sinva)) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',sa(1),',',sa(2),'],', &
                                                                       'inva=[',sinva(1),',',sinva(2),']')

         else

            !> Copy data in
            inva = a

            !> Compute matrix inverse
            call la_invert_c(inva,pivot=pivot,err=err0)

         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_invert_split_c

     ! Invert matrix in place
     function la_inverse_c(a,err) result(inva)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Output matrix inverse
         complex(sp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call la_invert_c(inva,err=err)

     end function la_inverse_c

     ! Inverse matrix operator
     function la_inverse_c_operator(a) result(inva)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Result matrix
         complex(sp),allocatable :: inva(:,:)

         type(la_state) :: err

         inva = la_inverse_c(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function la_inverse_c_operator

     ! Compute the in-place square matrix inverse of a
     subroutine la_invert_z(a,pivot,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout) :: a(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork,npiv
         integer(ilp),pointer :: ipiv(:)
         complex(dp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n .or. npiv < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = la_ilaenv(1,'zgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
         if (.not. present(pivot)) deallocate (ipiv)
         call err0%handle(err)

     end subroutine la_invert_z

     ! Compute the square matrix inverse of a into a second matrix
     subroutine la_invert_split_z(a,inva,pivot,err)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Inverse matrix inva[n,n]
         complex(dp),intent(out) :: inva(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: sa(2),sinva(2)
         character(*),parameter :: this = 'invert'

         sa = shape(a,kind=ilp)
         sinva = shape(inva,kind=ilp)

         if (any(sa /= sinva)) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',sa(1),',',sa(2),'],', &
                                                                       'inva=[',sinva(1),',',sinva(2),']')

         else

            !> Copy data in
            inva = a

            !> Compute matrix inverse
            call la_invert_z(inva,pivot=pivot,err=err0)

         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_invert_split_z

     ! Invert matrix in place
     function la_inverse_z(a,err) result(inva)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Output matrix inverse
         complex(dp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call la_invert_z(inva,err=err)

     end function la_inverse_z

     ! Inverse matrix operator
     function la_inverse_z_operator(a) result(inva)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Result matrix
         complex(dp),allocatable :: inva(:,:)

         type(la_state) :: err

         inva = la_inverse_z(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function la_inverse_z_operator

     ! Compute the in-place square matrix inverse of a
     subroutine la_invert_w(a,pivot,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout) :: a(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork,npiv
         integer(ilp),pointer :: ipiv(:)
         complex(qp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n .or. npiv < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = la_ilaenv(1,'wgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
         if (.not. present(pivot)) deallocate (ipiv)
         call err0%handle(err)

     end subroutine la_invert_w

     ! Compute the square matrix inverse of a into a second matrix
     subroutine la_invert_split_w(a,inva,pivot,err)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Inverse matrix inva[n,n]
         complex(qp),intent(out) :: inva(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: sa(2),sinva(2)
         character(*),parameter :: this = 'invert'

         sa = shape(a,kind=ilp)
         sinva = shape(inva,kind=ilp)

         if (any(sa /= sinva)) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',sa(1),',',sa(2),'],', &
                                                                       'inva=[',sinva(1),',',sinva(2),']')

         else

            !> Copy data in
            inva = a

            !> Compute matrix inverse
            call la_invert_w(inva,pivot=pivot,err=err0)

         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_invert_split_w

     ! Invert matrix in place
     function la_inverse_w(a,err) result(inva)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Output matrix inverse
         complex(qp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call la_invert_w(inva,err=err)

     end function la_inverse_w

     ! Inverse matrix operator
     function la_inverse_w_operator(a) result(inva)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Result matrix
         complex(qp),allocatable :: inva(:,:)

         type(la_state) :: err

         inva = la_inverse_w(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function la_inverse_w_operator

end module la_inverse
