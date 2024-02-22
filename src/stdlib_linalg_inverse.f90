
! Compute matrix inverse
module stdlib_linalg_inverse
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Function interface return the matrix inverse
     public :: inv
     !> Subroutine interface: invert matrix inplace
     public :: invert
     !> Operator interface: .inv.A returns the matrix inverse of A
     public :: operator(.inv.)

     ! Numpy: inv(a)
     ! Scipy: inv(a, overwrite_a=False, check_finite=True)
     ! IMSL: .i.a

     ! Function interface
     interface inv
        module procedure stdlib_linalg_inverse_s
        module procedure stdlib_linalg_inverse_d
        module procedure stdlib_linalg_inverse_q
        module procedure stdlib_linalg_inverse_c
        module procedure stdlib_linalg_inverse_z
        module procedure stdlib_linalg_inverse_w
     end interface inv

     ! Subroutine interface: in-place factorization
     interface invert
        module procedure stdlib_linalg_invert_s
        module procedure stdlib_linalg_invert_d
        module procedure stdlib_linalg_invert_q
        module procedure stdlib_linalg_invert_c
        module procedure stdlib_linalg_invert_z
        module procedure stdlib_linalg_invert_w
     end interface invert

     ! Operator interface
     interface operator(.inv.)
        module procedure stdlib_linalg_inverse_s_operator
        module procedure stdlib_linalg_inverse_d_operator
        module procedure stdlib_linalg_inverse_q_operator
        module procedure stdlib_linalg_inverse_c_operator
        module procedure stdlib_linalg_inverse_z_operator
        module procedure stdlib_linalg_inverse_w_operator
     end interface operator(.inv.)

     contains

     ! Compute the in-place square matrix inverse of a
     subroutine stdlib_linalg_invert_s(a,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         real(sp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            goto 1
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = stdlib_ilaenv(1,'sgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_invert_s

     ! Invert matrix in place
     function stdlib_linalg_inverse_s(a,err) result(inva)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Output matrix inverse
         real(sp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call stdlib_linalg_invert_s(inva,err)

     end function stdlib_linalg_inverse_s

     ! Inverse matrix operator
     function stdlib_linalg_inverse_s_operator(a) result(inva)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Result matrix
         real(sp),allocatable :: inva(:,:)

         type(linalg_state) :: err

         inva = stdlib_linalg_inverse_s(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function stdlib_linalg_inverse_s_operator

     ! Compute the in-place square matrix inverse of a
     subroutine stdlib_linalg_invert_d(a,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         real(dp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            goto 1
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = stdlib_ilaenv(1,'dgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_invert_d

     ! Invert matrix in place
     function stdlib_linalg_inverse_d(a,err) result(inva)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Output matrix inverse
         real(dp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call stdlib_linalg_invert_d(inva,err)

     end function stdlib_linalg_inverse_d

     ! Inverse matrix operator
     function stdlib_linalg_inverse_d_operator(a) result(inva)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Result matrix
         real(dp),allocatable :: inva(:,:)

         type(linalg_state) :: err

         inva = stdlib_linalg_inverse_d(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function stdlib_linalg_inverse_d_operator

     ! Compute the in-place square matrix inverse of a
     subroutine stdlib_linalg_invert_q(a,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         real(qp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            goto 1
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = stdlib_ilaenv(1,'qgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_invert_q

     ! Invert matrix in place
     function stdlib_linalg_inverse_q(a,err) result(inva)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Output matrix inverse
         real(qp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call stdlib_linalg_invert_q(inva,err)

     end function stdlib_linalg_inverse_q

     ! Inverse matrix operator
     function stdlib_linalg_inverse_q_operator(a) result(inva)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Result matrix
         real(qp),allocatable :: inva(:,:)

         type(linalg_state) :: err

         inva = stdlib_linalg_inverse_q(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function stdlib_linalg_inverse_q_operator

     ! Compute the in-place square matrix inverse of a
     subroutine stdlib_linalg_invert_c(a,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         complex(sp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            goto 1
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = stdlib_ilaenv(1,'cgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_invert_c

     ! Invert matrix in place
     function stdlib_linalg_inverse_c(a,err) result(inva)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Output matrix inverse
         complex(sp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call stdlib_linalg_invert_c(inva,err)

     end function stdlib_linalg_inverse_c

     ! Inverse matrix operator
     function stdlib_linalg_inverse_c_operator(a) result(inva)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Result matrix
         complex(sp),allocatable :: inva(:,:)

         type(linalg_state) :: err

         inva = stdlib_linalg_inverse_c(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function stdlib_linalg_inverse_c_operator

     ! Compute the in-place square matrix inverse of a
     subroutine stdlib_linalg_invert_z(a,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         complex(dp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            goto 1
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = stdlib_ilaenv(1,'zgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_invert_z

     ! Invert matrix in place
     function stdlib_linalg_inverse_z(a,err) result(inva)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Output matrix inverse
         complex(dp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call stdlib_linalg_invert_z(inva,err)

     end function stdlib_linalg_inverse_z

     ! Inverse matrix operator
     function stdlib_linalg_inverse_z_operator(a) result(inva)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Result matrix
         complex(dp),allocatable :: inva(:,:)

         type(linalg_state) :: err

         inva = stdlib_linalg_inverse_z(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function stdlib_linalg_inverse_z_operator

     ! Compute the in-place square matrix inverse of a
     subroutine stdlib_linalg_invert_w(a,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout) :: a(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,nb,lwork
         integer(ilp),allocatable :: ipiv(:)
         complex(qp),allocatable :: work(:)
         character(*),parameter :: this = 'invert'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (lda < 1 .or. n < 1 .or. lda /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=[',lda,',',n,']')
            goto 1
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Factorize matrix (overwrite result)
         call getrf(lda,n,a,lda,ipiv,info)

         ! Return codes from getrf and getri are identical
         if (info == 0) then

            ! Get optimal worksize (returned in work(1)) (apply 2% safety parameter)
            nb = stdlib_ilaenv(1,'wgetri',' ',n,-1,-1,-1)
            lwork = nint(1.02*n*nb,kind=ilp)

            allocate (work(lwork))

            ! Invert matrix
            call getri(n,a,lda,ipiv,work,lwork,info)

         end if

         select case (info)
            case (0)
                ! Success
            case (:-1)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                ! Matrix is singular
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output and return
1 call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_invert_w

     ! Invert matrix in place
     function stdlib_linalg_inverse_w(a,err) result(inva)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Output matrix inverse
         complex(qp),allocatable :: inva(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Allocate with copy
         allocate (inva,source=a)

         !> Compute matrix inverse
         call stdlib_linalg_invert_w(inva,err)

     end function stdlib_linalg_inverse_w

     ! Inverse matrix operator
     function stdlib_linalg_inverse_w_operator(a) result(inva)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Result matrix
         complex(qp),allocatable :: inva(:,:)

         type(linalg_state) :: err

         inva = stdlib_linalg_inverse_w(a,err)

         ! On error, return an empty matrix
         if (err%error()) then
            if (allocated(inva)) deallocate (inva)
            allocate (inva(0,0))
         end if

     end function stdlib_linalg_inverse_w_operator

end module stdlib_linalg_inverse
