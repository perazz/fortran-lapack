!> Determinant of a rectangular matrix
module la_determinant
     use la_constants
     use la_blas
     use la_lapack
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     public :: det

     character(*),parameter :: this = 'determinant'

     !> @brief Compute the determinant of a rectangular matrix.
     !!
     !! This function computes the determinant of a real or complex rectangular matrix \f$ A \f$.
     !! The determinant is calculated using LU factorization.
     !!
     !! @param[in,out] A The input square matrix of size \f$ [m, n] \f$. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] overwrite_a (Optional) If `.true.`, A may be overwritten and destroyed. Default is `.false.`.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided, 
     !!                 the function will stop execution.
     !!
     !! @return The determinant of the matrix \f$ A \f$. The result is a `real` scalar value of the same 
     !!         kind as the input matrix.
     !!
     !! @note This function relies on a matrix factorization approach (e.g., LU decomposition) to compute
     !!       the determinant efficiently from the [getrf](@ref la_lapack::getrf) backend.
     !! @warning If `overwrite_a` is enabled, the original contents of A will be lost.
     !!
     interface det
        module procedure la_sdeterminant
        module procedure la_ddeterminant
        module procedure la_qdeterminant
        module procedure la_cdeterminant
        module procedure la_zdeterminant
        module procedure la_wdeterminant
     end interface det

     contains

     ! Compute determinant of a square matrix A
     function la_sdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result: matrix determinant
         real(sp) :: det

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_sp
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         select case (m)
            case (0)
                ! Empty array has determinant 1 because math
                det = 1.0_sp

            case (1)
                ! Scalar
                det = a(1,1)

            case default

                ! Find determinant from LU decomposition

                ! Initialize a matrix temporary
                if (copy_a) then
                   allocate (amat(m,n),source=a)
                else
                   amat => a
                end if

                ! Pivot indices
                allocate (ipiv(n))

                ! Compute determinant from LU factorization, then calculate the product of
                ! all diagonal entries of the U factor.
                call getrf(m,n,amat,m,ipiv,info)

                select case (info)
                   case (0)
                       ! Success: compute determinant

                       ! Start with real 1.0
                       det = 1.0_sp
                       perm = 0
                       do k = 1,n
                          if (ipiv(k) /= k) perm = perm + 1
                          det = det*amat(k,k)
                       end do
                       if (mod(perm,2) /= 0) det = -det

                   case (:-1)
                       err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = la_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (copy_a) deallocate (amat)

         end select

         ! Process output and return
         call err0%handle(err)

     end function la_sdeterminant

     ! Compute determinant of a square matrix A
     function la_ddeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result: matrix determinant
         real(dp) :: det

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_dp
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         select case (m)
            case (0)
                ! Empty array has determinant 1 because math
                det = 1.0_dp

            case (1)
                ! Scalar
                det = a(1,1)

            case default

                ! Find determinant from LU decomposition

                ! Initialize a matrix temporary
                if (copy_a) then
                   allocate (amat(m,n),source=a)
                else
                   amat => a
                end if

                ! Pivot indices
                allocate (ipiv(n))

                ! Compute determinant from LU factorization, then calculate the product of
                ! all diagonal entries of the U factor.
                call getrf(m,n,amat,m,ipiv,info)

                select case (info)
                   case (0)
                       ! Success: compute determinant

                       ! Start with real 1.0
                       det = 1.0_dp
                       perm = 0
                       do k = 1,n
                          if (ipiv(k) /= k) perm = perm + 1
                          det = det*amat(k,k)
                       end do
                       if (mod(perm,2) /= 0) det = -det

                   case (:-1)
                       err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = la_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (copy_a) deallocate (amat)

         end select

         ! Process output and return
         call err0%handle(err)

     end function la_ddeterminant

     ! Compute determinant of a square matrix A
     function la_qdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result: matrix determinant
         real(qp) :: det

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_qp
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         select case (m)
            case (0)
                ! Empty array has determinant 1 because math
                det = 1.0_qp

            case (1)
                ! Scalar
                det = a(1,1)

            case default

                ! Find determinant from LU decomposition

                ! Initialize a matrix temporary
                if (copy_a) then
                   allocate (amat(m,n),source=a)
                else
                   amat => a
                end if

                ! Pivot indices
                allocate (ipiv(n))

                ! Compute determinant from LU factorization, then calculate the product of
                ! all diagonal entries of the U factor.
                call getrf(m,n,amat,m,ipiv,info)

                select case (info)
                   case (0)
                       ! Success: compute determinant

                       ! Start with real 1.0
                       det = 1.0_qp
                       perm = 0
                       do k = 1,n
                          if (ipiv(k) /= k) perm = perm + 1
                          det = det*amat(k,k)
                       end do
                       if (mod(perm,2) /= 0) det = -det

                   case (:-1)
                       err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = la_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (copy_a) deallocate (amat)

         end select

         ! Process output and return
         call err0%handle(err)

     end function la_qdeterminant

     ! Compute determinant of a square matrix A
     function la_cdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result: matrix determinant
         complex(sp) :: det

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_sp
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         select case (m)
            case (0)
                ! Empty array has determinant 1 because math
                det = 1.0_sp

            case (1)
                ! Scalar
                det = a(1,1)

            case default

                ! Find determinant from LU decomposition

                ! Initialize a matrix temporary
                if (copy_a) then
                   allocate (amat(m,n),source=a)
                else
                   amat => a
                end if

                ! Pivot indices
                allocate (ipiv(n))

                ! Compute determinant from LU factorization, then calculate the product of
                ! all diagonal entries of the U factor.
                call getrf(m,n,amat,m,ipiv,info)

                select case (info)
                   case (0)
                       ! Success: compute determinant

                       ! Start with real 1.0
                       det = 1.0_sp
                       perm = 0
                       do k = 1,n
                          if (ipiv(k) /= k) perm = perm + 1
                          det = det*amat(k,k)
                       end do
                       if (mod(perm,2) /= 0) det = -det

                   case (:-1)
                       err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = la_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (copy_a) deallocate (amat)

         end select

         ! Process output and return
         call err0%handle(err)

     end function la_cdeterminant

     ! Compute determinant of a square matrix A
     function la_zdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result: matrix determinant
         complex(dp) :: det

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_dp
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         select case (m)
            case (0)
                ! Empty array has determinant 1 because math
                det = 1.0_dp

            case (1)
                ! Scalar
                det = a(1,1)

            case default

                ! Find determinant from LU decomposition

                ! Initialize a matrix temporary
                if (copy_a) then
                   allocate (amat(m,n),source=a)
                else
                   amat => a
                end if

                ! Pivot indices
                allocate (ipiv(n))

                ! Compute determinant from LU factorization, then calculate the product of
                ! all diagonal entries of the U factor.
                call getrf(m,n,amat,m,ipiv,info)

                select case (info)
                   case (0)
                       ! Success: compute determinant

                       ! Start with real 1.0
                       det = 1.0_dp
                       perm = 0
                       do k = 1,n
                          if (ipiv(k) /= k) perm = perm + 1
                          det = det*amat(k,k)
                       end do
                       if (mod(perm,2) /= 0) det = -det

                   case (:-1)
                       err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = la_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (copy_a) deallocate (amat)

         end select

         ! Process output and return
         call err0%handle(err)

     end function la_zdeterminant

     ! Compute determinant of a square matrix A
     function la_wdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result: matrix determinant
         complex(qp) :: det

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_qp
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         select case (m)
            case (0)
                ! Empty array has determinant 1 because math
                det = 1.0_qp

            case (1)
                ! Scalar
                det = a(1,1)

            case default

                ! Find determinant from LU decomposition

                ! Initialize a matrix temporary
                if (copy_a) then
                   allocate (amat(m,n),source=a)
                else
                   amat => a
                end if

                ! Pivot indices
                allocate (ipiv(n))

                ! Compute determinant from LU factorization, then calculate the product of
                ! all diagonal entries of the U factor.
                call getrf(m,n,amat,m,ipiv,info)

                select case (info)
                   case (0)
                       ! Success: compute determinant

                       ! Start with real 1.0
                       det = 1.0_qp
                       perm = 0
                       do k = 1,n
                          if (ipiv(k) /= k) perm = perm + 1
                          det = det*amat(k,k)
                       end do
                       if (mod(perm,2) /= 0) det = -det

                   case (:-1)
                       err0 = la_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = la_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (copy_a) deallocate (amat)

         end select

         ! Process output and return
         call err0%handle(err)

     end function la_wdeterminant

end module la_determinant
