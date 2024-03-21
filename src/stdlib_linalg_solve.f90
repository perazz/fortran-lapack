module stdlib_linalg_solve
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Solve a linear system
     public :: solve

     ! Scipy: solve(a, b, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True, assume_a='gen', transposed=False)[source]#
     ! IMSL: lu_solve(a, b, transpose=False)

     interface solve
        module procedure stdlib_linalg_ssolve_one
        module procedure stdlib_linalg_dsolve_one
        module procedure stdlib_linalg_qsolve_one
        module procedure stdlib_linalg_csolve_one
        module procedure stdlib_linalg_zsolve_one
        module procedure stdlib_linalg_wsolve_one
        module procedure stdlib_linalg_ssolve_multiple
        module procedure stdlib_linalg_dsolve_multiple
        module procedure stdlib_linalg_qsolve_multiple
        module procedure stdlib_linalg_csolve_multiple
        module procedure stdlib_linalg_zsolve_multiple
        module procedure stdlib_linalg_wsolve_multiple
     end interface solve

     contains

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_ssolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_ssolve_one

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_dsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_dsolve_one

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_qsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_qsolve_one

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_csolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_csolve_one

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_zsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_zsolve_one

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_wsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_wsolve_one

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_ssolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_ssolve_multiple

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_dsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_dsolve_multiple

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_qsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_qsolve_multiple

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_csolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_csolve_multiple

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_zsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_zsolve_multiple

     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_wsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'solve'

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            goto 1
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Pivot indices
         allocate (ipiv(n))

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         allocate (x,source=b)
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         if (.not. copy_a) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_wsolve_multiple

end module stdlib_linalg_solve
