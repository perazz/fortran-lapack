
module stdlib_linalg_determinant
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Determinant of a rectangular matrix
     public :: det

     character(*),parameter :: this = 'determinant'

     ! Numpy: det(a)
     ! Scipy: det(a, overwrite_a=False, check_finite=True)
     ! IMSL: DET(a)

     interface det
        module procedure stdlib_linalg_sdeterminant
        module procedure stdlib_linalg_ddeterminant
        module procedure stdlib_linalg_qdeterminant
        module procedure stdlib_linalg_cdeterminant
        module procedure stdlib_linalg_zdeterminant
        module procedure stdlib_linalg_wdeterminant
     end interface det

     contains

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_sdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result: matrix determinant
         real(sp) :: det

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_sp
            goto 1
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
                       err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (.not. copy_a) deallocate (amat)

         end select

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_sdeterminant

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_ddeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result: matrix determinant
         real(dp) :: det

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_dp
            goto 1
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
                       err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (.not. copy_a) deallocate (amat)

         end select

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_ddeterminant

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_qdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result: matrix determinant
         real(qp) :: det

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_qp
            goto 1
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
                       err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (.not. copy_a) deallocate (amat)

         end select

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_qdeterminant

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_cdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result: matrix determinant
         complex(sp) :: det

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_sp
            goto 1
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
                       err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (.not. copy_a) deallocate (amat)

         end select

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_cdeterminant

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_zdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result: matrix determinant
         complex(dp) :: det

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_dp
            goto 1
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
                       err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (.not. copy_a) deallocate (amat)

         end select

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_zdeterminant

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_wdeterminant(a,overwrite_a,err) result(det)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result: matrix determinant
         complex(qp) :: det

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: m,n,info,perm,k
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: amat(:,:)

         !> Matrix determinant size
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m /= n .or. .not. min(m,n) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid or non-square matrix: a=[',m,',',n,']')
            det = 0.0_qp
            goto 1
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
                       err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',m,',',n,']')
                   case (1:)
                       err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
                   case default
                       err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
                end select

                if (.not. copy_a) deallocate (amat)

         end select

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_wdeterminant

end module stdlib_linalg_determinant
