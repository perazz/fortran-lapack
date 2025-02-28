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


     !> @brief Compute the inverse of a square matrix in-place.
     !!
     !! This subroutine computes the inverse of a real or complex square matrix \f$ A \f$ in-place.
     !! The inverse is computed using an LU decomposition with partial pivoting.
     !!
     !! @param[in,out] A The input square matrix of size \f$ [n, n] \f$. It is replaced by its inverse \f$ A^{-1} \f$.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @note This subroutine is useful when memory efficiency is a priority, as it avoids additional allocations.
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

     ! Subroutine interface: in-place factorization
     interface invert
        module procedure la_invert_s
        module procedure la_invert_d
        module procedure la_invert_q
        module procedure la_invert_c
        module procedure la_invert_z
        module procedure la_invert_w
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
     subroutine la_invert_s(a,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         real(sp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            call err0%handle(err)
            return
         end if

         ! Pivot indices
         allocate (ipiv(n))

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
         call err0%handle(err)

     end subroutine la_invert_s

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
         call la_invert_s(inva,err)

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
     subroutine la_invert_d(a,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         real(dp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            call err0%handle(err)
            return
         end if

         ! Pivot indices
         allocate (ipiv(n))

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
         call err0%handle(err)

     end subroutine la_invert_d

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
         call la_invert_d(inva,err)

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
     subroutine la_invert_q(a,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         real(qp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            call err0%handle(err)
            return
         end if

         ! Pivot indices
         allocate (ipiv(n))

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
         call err0%handle(err)

     end subroutine la_invert_q

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
         call la_invert_q(inva,err)

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
     subroutine la_invert_c(a,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         complex(sp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            call err0%handle(err)
            return
         end if

         ! Pivot indices
         allocate (ipiv(n))

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
         call err0%handle(err)

     end subroutine la_invert_c

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
         call la_invert_c(inva,err)

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
     subroutine la_invert_z(a,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         complex(dp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            call err0%handle(err)
            return
         end if

         ! Pivot indices
         allocate (ipiv(n))

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
         call err0%handle(err)

     end subroutine la_invert_z

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
         call la_invert_z(inva,err)

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
     subroutine la_invert_w(a,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         complex(qp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            call err0%handle(err)
            return
         end if

         ! Pivot indices
         allocate (ipiv(n))

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
         call err0%handle(err)

     end subroutine la_invert_w

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
         call la_invert_w(inva,err)

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
