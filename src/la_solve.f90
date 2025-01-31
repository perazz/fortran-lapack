module la_solve
     use la_constants
     use la_blas
     use la_lapack
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> @brief Solve a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !! \f$ A \cdot X = B \f$
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     public :: solve

     interface solve
        module procedure la_ssolve_one
        module procedure la_dsolve_one
        module procedure la_qsolve_one
        module procedure la_csolve_one
        module procedure la_zsolve_one
        module procedure la_wsolve_one
        module procedure la_ssolve_multiple
        module procedure la_dsolve_multiple
        module procedure la_qsolve_multiple
        module procedure la_csolve_multiple
        module procedure la_zsolve_multiple
        module procedure la_wsolve_multiple
     end interface solve
     
     character(*),parameter :: this = 'solve'

     contains
     
     elemental subroutine handle_gesv_info(info,lda,n,nrhs,err)
         integer(ilp),intent(in) :: info,lda,n,nrhs
         type(la_state),intent(out) :: err

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [lda,n])
            case (-7)
                err = la_state(this,LINALG_ERROR,'invalid matrix size a=', [lda,n])
            case (1:)
                err = la_state(this,LINALG_ERROR,'singular matrix')
            case default
                err = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gesv_info

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_ssolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_ssolve_one

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_dsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_dsolve_one

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_qsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_qsolve_one

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_csolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_csolve_one

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_zsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_zsolve_one

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_wsolve_one(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_wsolve_one

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_ssolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_ssolve_multiple

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_dsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_dsolve_multiple

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_qsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_qsolve_multiple

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_csolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_csolve_multiple

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_zsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_zsolve_multiple

     !> @brief Solves a system of linear equations A * X = B.
     !!
     !! This function computes the solution to a real system of linear equations:
     !!
     !!     A * X = B
     !!
     !! where A is an `n x n` square matrix, and B is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution X is returned as an allocatable array.
     !!
     !! @tparam T Data type of the matrix and vectors (real or complex).
     !! @tparam N Size of the square matrix A.
     !! @tparam NRHS Number of right-hand sides (optional, default is 1).
     !!
     !! @param[in,out] A The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] B The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return Solution matrix X of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK LU decomposition based solvers `*GESV`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     !! @example
     !!
     !!     real(kind=dp), allocatable :: A(:,:), B(:), X(:)
     !!     allocate(A(3,3), B(3))
     !!     A = reshape([3.0_dp, 2.0_dp, -1.0_dp, &
     !!                  2.0_dp, -2.0_dp, 4.0_dp, &
     !!                 -1.0_dp, 0.5_dp, -1.0_dp], [3,3])
     !!     B = [1.0_dp, -2.0_dp, 0.0_dp]
     !!     X = solve(A, B)
     !!
     function la_wsolve_multiple(a,b,overwrite_a,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp),allocatable :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. lda /= n .or. ldb /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,']')
            allocate (x(0,0))
            call err0%handle(err)
            return
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
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end function la_wsolve_multiple

end module la_solve
