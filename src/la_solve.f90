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

     !> @brief Solve a system of linear equations A * X = B into a pre-allocated array.
     !!
     !! ### Summary
     !! Subroutine interface for solving a linear system through an LU factorization.
     !!
     !! ### Description
     !!
     !! This interface solves the linear system
     !!
     !! \f$ A \cdot X = B \f$
     !!
     !! writing the result into the caller's array `x` instead of returning a new one.
     !! Storage for the pivot indices may also be provided; when both `x` and `pivot` are
     !! given and `overwrite_a` is set, the call performs no internal allocation.
     !! One (`b(:)`) or many (`b(:,:)`) right-hand sides are solved at once.
     !!
     !! @note The solution is based on LAPACK's LU decomposition solvers [GESV](@ref la_lapack::gesv).
     !!
     !! @param[in,out] a The input square matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] b The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in,out] x The solution vector (size `n`) or matrix (size `n x nrhs`), overwritten on return.
     !! @param[in,out] pivot (Optional) Storage array for the `n` diagonal pivot indices.
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the subroutine will stop execution.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     public :: solve_lu

     !> @brief Solve a Hermitian positive definite system A * X = B.
     !!
     !! ### Summary
     !! Subroutine interface for factorizing and solving a Cholesky system in one call.
     !!
     !! ### Description
     !!
     !! This interface computes the Cholesky factorization of the `real` symmetric or `complex`
     !! Hermitian positive definite matrix A and solves
     !!
     !! \f$ A \cdot X = B \f$
     !!
     !! for one (`b(:)`) or many (`b(:,:)`) right-hand sides, writing the result into the caller's
     !! array `x`. Only the triangle `lower` selects is read.
     !!
     !! @note The solution is based on LAPACK's [POSV](@ref la_lapack::posv) drivers.
     !!
     !! @param[in,out] a The input matrix of size `n x n`. If `overwrite_a` is true, it is
     !!                  overwritten with its Cholesky factor.
     !! @param[in] b The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in,out] x The solution vector (size `n`) or matrix (size `n x nrhs`), overwritten on return.
     !! @param[in] lower (Optional) If true, the lower triangle of A is used. Default is true.
     !! @param[in] overwrite_a (Optional) If true, A may be overwritten and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the subroutine will stop execution.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A may be lost.
     !!
     public :: solve_chol

     !> @brief Solve A * X = B from a pre-computed lower Cholesky factor.
     !!
     !! ### Summary
     !! Subroutine interface for solving with the lower triangular factor of \f$ A = L \cdot L^H \f$.
     !!
     !! ### Description
     !!
     !! This interface solves \f$ A \cdot X = B \f$ for one or many right-hand sides, given the
     !! lower Cholesky factor of A as returned by [cholesky](@ref la_cholesky::cholesky) with
     !! `lower=.true.`. The factor is not recomputed, so a repeated solve of the same system costs
     !! two triangular solves.
     !!
     !! @note The solution is based on LAPACK's [POTRS](@ref la_lapack::potrs) routines.
     !!
     !! @param[in] l The lower Cholesky factor of size `n x n`.
     !! @param[in] b The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in,out] x The solution vector (size `n`) or matrix (size `n x nrhs`), overwritten on return.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the subroutine will stop execution.
     !!
     public :: solve_lower_chol

     !> @brief Solve A * X = B from a pre-computed upper Cholesky factor.
     !!
     !! ### Summary
     !! Subroutine interface for solving with the upper triangular factor of \f$ A = U^H \cdot U \f$.
     !!
     !! ### Description
     !!
     !! This interface solves \f$ A \cdot X = B \f$ for one or many right-hand sides, given the
     !! upper Cholesky factor of A as returned by [cholesky](@ref la_cholesky::cholesky) with
     !! `lower=.false.`. The factor is not recomputed, so a repeated solve of the same system costs
     !! two triangular solves.
     !!
     !! @note The solution is based on LAPACK's [POTRS](@ref la_lapack::potrs) routines.
     !!
     !! @param[in] u The upper Cholesky factor of size `n x n`.
     !! @param[in] b The right-hand side vector (size `n`) or matrix (size `n x nrhs`).
     !! @param[in,out] x The solution vector (size `n`) or matrix (size `n x nrhs`), overwritten on return.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the subroutine will stop execution.
     !!
     public :: solve_upper_chol

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

     interface solve_lu
        module procedure la_ssolve_lu_one
        module procedure la_dsolve_lu_one
        module procedure la_qsolve_lu_one
        module procedure la_csolve_lu_one
        module procedure la_zsolve_lu_one
        module procedure la_wsolve_lu_one
        module procedure la_ssolve_lu_multiple
        module procedure la_dsolve_lu_multiple
        module procedure la_qsolve_lu_multiple
        module procedure la_csolve_lu_multiple
        module procedure la_zsolve_lu_multiple
        module procedure la_wsolve_lu_multiple
     end interface solve_lu

     interface solve_chol
        module procedure la_ssolve_chol_one
        module procedure la_dsolve_chol_one
        module procedure la_qsolve_chol_one
        module procedure la_csolve_chol_one
        module procedure la_zsolve_chol_one
        module procedure la_wsolve_chol_one
        module procedure la_ssolve_chol_multiple
        module procedure la_dsolve_chol_multiple
        module procedure la_qsolve_chol_multiple
        module procedure la_csolve_chol_multiple
        module procedure la_zsolve_chol_multiple
        module procedure la_wsolve_chol_multiple
     end interface solve_chol

     interface solve_lower_chol
        module procedure la_ssolve_lower_chol_one
        module procedure la_dsolve_lower_chol_one
        module procedure la_qsolve_lower_chol_one
        module procedure la_csolve_lower_chol_one
        module procedure la_zsolve_lower_chol_one
        module procedure la_wsolve_lower_chol_one
        module procedure la_ssolve_lower_chol_multiple
        module procedure la_dsolve_lower_chol_multiple
        module procedure la_qsolve_lower_chol_multiple
        module procedure la_csolve_lower_chol_multiple
        module procedure la_zsolve_lower_chol_multiple
        module procedure la_wsolve_lower_chol_multiple
     end interface solve_lower_chol

     interface solve_upper_chol
        module procedure la_ssolve_upper_chol_one
        module procedure la_dsolve_upper_chol_one
        module procedure la_qsolve_upper_chol_one
        module procedure la_csolve_upper_chol_one
        module procedure la_zsolve_upper_chol_one
        module procedure la_wsolve_upper_chol_one
        module procedure la_ssolve_upper_chol_multiple
        module procedure la_dsolve_upper_chol_multiple
        module procedure la_qsolve_upper_chol_multiple
        module procedure la_csolve_upper_chol_multiple
        module procedure la_zsolve_upper_chol_multiple
        module procedure la_wsolve_upper_chol_multiple
     end interface solve_upper_chol
     
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

     elemental subroutine handle_chol_info(info,triangle,n,nrhs,lda,ldb,err)
         integer(ilp),intent(in) :: info,n,nrhs,lda,ldb
         character,intent(in) :: triangle
         type(la_state),intent(out) :: err

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(this,LINALG_INTERNAL_ERROR,'invalid triangle selection: ', &
                                    triangle,'. should be U/L')
            case (-2)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-3)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid rhs size nrhs=',nrhs)
            case (-5)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid lda=',lda,': should be >=',n)
            case (-7)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid ldb=',ldb,': should be >=',n)
            case (1:)
                err = la_state(this,LINALG_ERROR,'matrix is not positive definite: ', &
                                    'leading minor of order',info,'is not positive definite')
            case default
                err = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_chol_info

     !> Linear system solve, _one, real(sp)
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

     !> Linear system solve into a pre-allocated array, _one, real(sp)
     pure subroutine la_ssolve_lu_one(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_lu_one

     !> Cholesky factorize and solve in one call, _one, real(sp)
     pure subroutine la_ssolve_chol_one(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         real(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_chol_one

     !> Solve from a pre-computed Cholesky factor, _one, real(sp)
     pure subroutine la_schol_solve_one(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         real(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         real(sp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_schol_solve_one

     !> Solve from a pre-computed lower Cholesky factor, _one, real(sp)
     pure subroutine la_ssolve_lower_chol_one(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         real(sp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_schol_solve_one(l,b,x,'L',err)

     end subroutine la_ssolve_lower_chol_one

     !> Solve from a pre-computed upper Cholesky factor, _one, real(sp)
     pure subroutine la_ssolve_upper_chol_one(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         real(sp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_schol_solve_one(u,b,x,'U',err)

     end subroutine la_ssolve_upper_chol_one

     !> Linear system solve, _one, real(dp)
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

     !> Linear system solve into a pre-allocated array, _one, real(dp)
     pure subroutine la_dsolve_lu_one(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_lu_one

     !> Cholesky factorize and solve in one call, _one, real(dp)
     pure subroutine la_dsolve_chol_one(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         real(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_chol_one

     !> Solve from a pre-computed Cholesky factor, _one, real(dp)
     pure subroutine la_dchol_solve_one(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         real(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         real(dp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dchol_solve_one

     !> Solve from a pre-computed lower Cholesky factor, _one, real(dp)
     pure subroutine la_dsolve_lower_chol_one(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         real(dp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_dchol_solve_one(l,b,x,'L',err)

     end subroutine la_dsolve_lower_chol_one

     !> Solve from a pre-computed upper Cholesky factor, _one, real(dp)
     pure subroutine la_dsolve_upper_chol_one(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         real(dp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_dchol_solve_one(u,b,x,'U',err)

     end subroutine la_dsolve_upper_chol_one

     !> Linear system solve, _one, real(qp)
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

     !> Linear system solve into a pre-allocated array, _one, real(qp)
     pure subroutine la_qsolve_lu_one(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_lu_one

     !> Cholesky factorize and solve in one call, _one, real(qp)
     pure subroutine la_qsolve_chol_one(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         real(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_chol_one

     !> Solve from a pre-computed Cholesky factor, _one, real(qp)
     pure subroutine la_qchol_solve_one(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         real(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         real(qp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qchol_solve_one

     !> Solve from a pre-computed lower Cholesky factor, _one, real(qp)
     pure subroutine la_qsolve_lower_chol_one(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         real(qp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_qchol_solve_one(l,b,x,'L',err)

     end subroutine la_qsolve_lower_chol_one

     !> Solve from a pre-computed upper Cholesky factor, _one, real(qp)
     pure subroutine la_qsolve_upper_chol_one(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         real(qp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_qchol_solve_one(u,b,x,'U',err)

     end subroutine la_qsolve_upper_chol_one

     !> Linear system solve, _one, complex(sp)
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

     !> Linear system solve into a pre-allocated array, _one, complex(sp)
     pure subroutine la_csolve_lu_one(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_lu_one

     !> Cholesky factorize and solve in one call, _one, complex(sp)
     pure subroutine la_csolve_chol_one(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         complex(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_chol_one

     !> Solve from a pre-computed Cholesky factor, _one, complex(sp)
     pure subroutine la_cchol_solve_one(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         complex(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         complex(sp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_cchol_solve_one

     !> Solve from a pre-computed lower Cholesky factor, _one, complex(sp)
     pure subroutine la_csolve_lower_chol_one(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         complex(sp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_cchol_solve_one(l,b,x,'L',err)

     end subroutine la_csolve_lower_chol_one

     !> Solve from a pre-computed upper Cholesky factor, _one, complex(sp)
     pure subroutine la_csolve_upper_chol_one(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         complex(sp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_cchol_solve_one(u,b,x,'U',err)

     end subroutine la_csolve_upper_chol_one

     !> Linear system solve, _one, complex(dp)
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

     !> Linear system solve into a pre-allocated array, _one, complex(dp)
     pure subroutine la_zsolve_lu_one(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_lu_one

     !> Cholesky factorize and solve in one call, _one, complex(dp)
     pure subroutine la_zsolve_chol_one(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         complex(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_chol_one

     !> Solve from a pre-computed Cholesky factor, _one, complex(dp)
     pure subroutine la_zchol_solve_one(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         complex(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         complex(dp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zchol_solve_one

     !> Solve from a pre-computed lower Cholesky factor, _one, complex(dp)
     pure subroutine la_zsolve_lower_chol_one(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         complex(dp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_zchol_solve_one(l,b,x,'L',err)

     end subroutine la_zsolve_lower_chol_one

     !> Solve from a pre-computed upper Cholesky factor, _one, complex(dp)
     pure subroutine la_zsolve_upper_chol_one(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         complex(dp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_zchol_solve_one(u,b,x,'U',err)

     end subroutine la_zsolve_upper_chol_one

     !> Linear system solve, _one, complex(qp)
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

     !> Linear system solve into a pre-allocated array, _one, complex(qp)
     pure subroutine la_wsolve_lu_one(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_lu_one

     !> Cholesky factorize and solve in one call, _one, complex(qp)
     pure subroutine la_wsolve_chol_one(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         complex(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_chol_one

     !> Solve from a pre-computed Cholesky factor, _one, complex(qp)
     pure subroutine la_wchol_solve_one(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         complex(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         complex(qp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wchol_solve_one

     !> Solve from a pre-computed lower Cholesky factor, _one, complex(qp)
     pure subroutine la_wsolve_lower_chol_one(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         complex(qp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_wchol_solve_one(l,b,x,'L',err)

     end subroutine la_wsolve_lower_chol_one

     !> Solve from a pre-computed upper Cholesky factor, _one, complex(qp)
     pure subroutine la_wsolve_upper_chol_one(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         complex(qp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_wchol_solve_one(u,b,x,'U',err)

     end subroutine la_wsolve_upper_chol_one

     !> Linear system solve, _multiple, real(sp)
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

     !> Linear system solve into a pre-allocated array, _multiple, real(sp)
     pure subroutine la_ssolve_lu_multiple(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         real(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_lu_multiple

     !> Cholesky factorize and solve in one call, _multiple, real(sp)
     pure subroutine la_ssolve_chol_multiple(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         real(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_chol_multiple

     !> Solve from a pre-computed Cholesky factor, _multiple, real(sp)
     pure subroutine la_schol_solve_multiple(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         real(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:,:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         real(sp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_schol_solve_multiple

     !> Solve from a pre-computed lower Cholesky factor, _multiple, real(sp)
     pure subroutine la_ssolve_lower_chol_multiple(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         real(sp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_schol_solve_multiple(l,b,x,'L',err)

     end subroutine la_ssolve_lower_chol_multiple

     !> Solve from a pre-computed upper Cholesky factor, _multiple, real(sp)
     pure subroutine la_ssolve_upper_chol_multiple(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         real(sp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_schol_solve_multiple(u,b,x,'U',err)

     end subroutine la_ssolve_upper_chol_multiple

     !> Linear system solve, _multiple, real(dp)
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

     !> Linear system solve into a pre-allocated array, _multiple, real(dp)
     pure subroutine la_dsolve_lu_multiple(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         real(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_lu_multiple

     !> Cholesky factorize and solve in one call, _multiple, real(dp)
     pure subroutine la_dsolve_chol_multiple(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         real(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_chol_multiple

     !> Solve from a pre-computed Cholesky factor, _multiple, real(dp)
     pure subroutine la_dchol_solve_multiple(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         real(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:,:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         real(dp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dchol_solve_multiple

     !> Solve from a pre-computed lower Cholesky factor, _multiple, real(dp)
     pure subroutine la_dsolve_lower_chol_multiple(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         real(dp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_dchol_solve_multiple(l,b,x,'L',err)

     end subroutine la_dsolve_lower_chol_multiple

     !> Solve from a pre-computed upper Cholesky factor, _multiple, real(dp)
     pure subroutine la_dsolve_upper_chol_multiple(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         real(dp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_dchol_solve_multiple(u,b,x,'U',err)

     end subroutine la_dsolve_upper_chol_multiple

     !> Linear system solve, _multiple, real(qp)
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

     !> Linear system solve into a pre-allocated array, _multiple, real(qp)
     pure subroutine la_qsolve_lu_multiple(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         real(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_lu_multiple

     !> Cholesky factorize and solve in one call, _multiple, real(qp)
     pure subroutine la_qsolve_chol_multiple(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         real(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_chol_multiple

     !> Solve from a pre-computed Cholesky factor, _multiple, real(qp)
     pure subroutine la_qchol_solve_multiple(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         real(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:,:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         real(qp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qchol_solve_multiple

     !> Solve from a pre-computed lower Cholesky factor, _multiple, real(qp)
     pure subroutine la_qsolve_lower_chol_multiple(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         real(qp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_qchol_solve_multiple(l,b,x,'L',err)

     end subroutine la_qsolve_lower_chol_multiple

     !> Solve from a pre-computed upper Cholesky factor, _multiple, real(qp)
     pure subroutine la_qsolve_upper_chol_multiple(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         real(qp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_qchol_solve_multiple(u,b,x,'U',err)

     end subroutine la_qsolve_upper_chol_multiple

     !> Linear system solve, _multiple, complex(sp)
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

     !> Linear system solve into a pre-allocated array, _multiple, complex(sp)
     pure subroutine la_csolve_lu_multiple(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         complex(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_lu_multiple

     !> Cholesky factorize and solve in one call, _multiple, complex(sp)
     pure subroutine la_csolve_chol_multiple(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         complex(sp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_chol_multiple

     !> Solve from a pre-computed Cholesky factor, _multiple, complex(sp)
     pure subroutine la_cchol_solve_multiple(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         complex(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:,:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         complex(sp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_cchol_solve_multiple

     !> Solve from a pre-computed lower Cholesky factor, _multiple, complex(sp)
     pure subroutine la_csolve_lower_chol_multiple(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         complex(sp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_cchol_solve_multiple(l,b,x,'L',err)

     end subroutine la_csolve_lower_chol_multiple

     !> Solve from a pre-computed upper Cholesky factor, _multiple, complex(sp)
     pure subroutine la_csolve_upper_chol_multiple(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         complex(sp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_cchol_solve_multiple(u,b,x,'U',err)

     end subroutine la_csolve_upper_chol_multiple

     !> Linear system solve, _multiple, complex(dp)
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

     !> Linear system solve into a pre-allocated array, _multiple, complex(dp)
     pure subroutine la_zsolve_lu_multiple(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         complex(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_lu_multiple

     !> Cholesky factorize and solve in one call, _multiple, complex(dp)
     pure subroutine la_zsolve_chol_multiple(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         complex(dp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_chol_multiple

     !> Solve from a pre-computed Cholesky factor, _multiple, complex(dp)
     pure subroutine la_zchol_solve_multiple(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         complex(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:,:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         complex(dp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zchol_solve_multiple

     !> Solve from a pre-computed lower Cholesky factor, _multiple, complex(dp)
     pure subroutine la_zsolve_lower_chol_multiple(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         complex(dp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_zchol_solve_multiple(l,b,x,'L',err)

     end subroutine la_zsolve_lower_chol_multiple

     !> Solve from a pre-computed upper Cholesky factor, _multiple, complex(dp)
     pure subroutine la_zsolve_upper_chol_multiple(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         complex(dp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_zchol_solve_multiple(u,b,x,'U',err)

     end subroutine la_zsolve_upper_chol_multiple

     !> Linear system solve, _multiple, complex(qp)
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

     !> Linear system solve into a pre-allocated array, _multiple, complex(qp)
     pure subroutine la_wsolve_lu_multiple(a,b,x,pivot,overwrite_a,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Storage array for the diagonal pivot indices
         integer(ilp),optional,intent(inout),target :: pivot(:)
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,npiv,info
         integer(ilp),pointer :: ipiv(:)
         logical(lk) :: copy_a
         complex(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Has a pre-allocated pivot storage array been provided?
         if (present(pivot)) then
            ipiv => pivot
         else
            allocate (ipiv(n))
         end if
         npiv = size(ipiv,kind=ilp)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs .or. npiv /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,'],', &
                                                                       'pivot=',npiv)
            if (.not. present(pivot)) deallocate (ipiv)
            call err0%handle(err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve system
         call gesv(n,nrhs,amat,lda,ipiv,xmat,ldb,info)

         ! Process output
         call handle_gesv_info(info,lda,n,nrhs,err0)

         if (copy_a) deallocate (amat)
         if (.not. present(pivot)) deallocate (ipiv)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_lu_multiple

     !> Cholesky factorize and solve in one call, _multiple, complex(qp)
     pure subroutine la_wsolve_chol_multiple(a,b,x,lower,overwrite_a,err)
         !> Input Hermitian positive definite matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] Use the lower triangular factorization? Default: .true.
         logical(lk),optional,intent(in) :: lower
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         logical(lk) :: lower_,copy_a
         character :: uplo
         complex(qp),pointer :: xmat(:,:),amat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         ! Default: use the lower triangle
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true._lk
         end if
         uplo = merge('L','U',lower_)

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize a matrix temporary
         if (copy_a) then
            allocate (amat(lda,n),source=a)
         else
            amat => a
         end if

         ! Initialize solution with the rhs: posv overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Factorize and solve
         call posv(uplo,n,nrhs,amat,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         if (copy_a) deallocate (amat)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_chol_multiple

     !> Solve from a pre-computed Cholesky factor, _multiple, complex(qp)
     pure subroutine la_wchol_solve_multiple(a,b,x,uplo,err)
         !> Cholesky factor a[n,n], lower or upper as uplo selects
         complex(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:,:)
         !> Triangle selector: 'L' for lower, 'U' for upper
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldb,ldx,nrhs,nrhsx,info
         complex(qp),pointer :: xmat(:,:)

         !> Problem sizes
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/max(1_ilp,ldb)
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/max(1_ilp,ldx)

         if (any([lda,n,ldb] < 1) .or. any([lda,ldb,ldx] /= n) .or. nrhsx /= nrhs) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'],', &
                                                                       'b=[',ldb,',',nrhs,'],', &
                                                                       'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            return
         end if

         ! Initialize solution with the rhs: potrs overwrites it with the solution
         x = b
         xmat(1:n,1:nrhs) => x

         ! Solve with the triangular factors
         call potrs(uplo,n,nrhs,a,lda,xmat,n,info)

         ! Process output
         call handle_chol_info(info,uplo,n,nrhs,lda,n,err0)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wchol_solve_multiple

     !> Solve from a pre-computed lower Cholesky factor, _multiple, complex(qp)
     pure subroutine la_wsolve_lower_chol_multiple(l,b,x,err)
         !> Lower Cholesky factor l[n,n] from cholesky(...,lower=.true.)
         complex(qp),intent(in) :: l(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_wchol_solve_multiple(l,b,x,'L',err)

     end subroutine la_wsolve_lower_chol_multiple

     !> Solve from a pre-computed upper Cholesky factor, _multiple, complex(qp)
     pure subroutine la_wsolve_upper_chol_multiple(u,b,x,err)
         !> Upper Cholesky factor u[n,n] from cholesky(...,lower=.false.)
         complex(qp),intent(in) :: u(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         call la_wchol_solve_multiple(u,b,x,'U',err)

     end subroutine la_wsolve_upper_chol_multiple

end module la_solve
