!> Least squares solution interface
module la_least_squares
     use la_constants
     use la_blas
     use la_lapack
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> @brief Compute a least squares solution to system \f$ A \cdot x = b \f$,
     !!  i.e. such that the 2-norm \f$ \|b - A \cdot x\| \f$ is minimized.
     !!
     !! This function computes the least-squares solution to a real system of linear equations:
     !!
     !! \f$ A \cdot x = b \f$
     !!
     !! where A is an `n x n` matrix and b is either a vector (`n`) or a matrix (`n x nrhs`).
     !! The solution `x` is returned as an allocatable array.
     !!
     !! @param[in,out] a The input matrix of size `n x n`. If `overwrite_a` is true,
     !!                  the contents of A may be modified during computation.
     !! @param[in] b The right-hand side vector (`n`) or matrix (`n x nrhs`).
     !! @param[in] cond (Optional) A cutoff for rank evaluation: singular values \f$ s(i) \f$ such that
     !!                \f$ s(i) \leq \text{cond} \cdot \max(s) \f$ are considered zero.
     !! @param[in] overwrite_a (Optional) If true, A and B may be overwritten and destroyed. Default is false.
     !! @param[out] rank (Optional) The rank of the matrix A.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     !! @return Solution matrix `x` of size `n` (for a single right-hand side) or `n x nrhs`.
     !!
     !! @note This function relies on LAPACK least-squares solvers such as `[*GELSS](@ref la_lapack::gelss)`.
     !!
     !! @warning If `overwrite_a` is enabled, the original contents of A and B may be lost.
     !!
     public :: lstsq

     !> @brief Compute the working space needed by a least squares solution of \f$ A \cdot x = b \f$.
     !!
     !! This subroutine returns the sizes of the real, integer and (for complex data) complex
     !! working arrays that [`solve_lstsq`](@ref la_least_squares::solve_lstsq) needs for the given
     !! problem, so that a repeated solve of problems of the same size performs no allocation.
     !!
     !! @param[in] a The input matrix of size \f$ [m,n] \f$. Only its shape is used.
     !! @param[in] b The right-hand side vector (\f$m\f$) or matrix (\f$m \times nrhs\f$). Only its shape is used.
     !! @param[out] lrwork The size of the real working array.
     !! @param[out] liwork The size of the integer working array.
     !! @param[out] lcwork The size of the complex working array. Present for complex data only.
     !!
     public :: lstsq_space

     !> @brief Compute a least squares solution of \f$ A \cdot x = b \f$ into a pre-allocated array.
     !!
     !! This subroutine computes the least-squares solution of a linear matrix problem and writes it
     !! into the caller's array, so that neither the solution nor the working space needs to be
     !! allocated internally.
     !!
     !! @param[in,out] a The input matrix of size \f$ [m,n] \f$. If `overwrite_a` is true,
     !!                  the contents of `a` may be modified during computation.
     !! @param[in] b The right-hand side vector (\f$m\f$) or matrix (\f$m \times nrhs\f$).
     !! @param[in,out] x The solution vector (\f$\ge n\f$) or matrix (\f$\ge n \times nrhs\f$).
     !! @param[in,out] real_storage (Optional) Pre-allocated real working space. Its size is checked with [`lstsq_space`](@ref la_least_squares::lstsq_space).
     !! @param[in,out] int_storage (Optional) Pre-allocated integer working space.
     !! @param[in,out] cmpl_storage (Optional) Pre-allocated complex working space. Complex data only.
     !! @param[in] cond (Optional) A cutoff for rank evaluation: singular values \f$ s(i) \f$ such that
     !!                \f$ s(i) \leq \text{cond} \cdot \max(s) \f$ are considered zero.
     !! @param[out] singvals (Optional) The \f$ \min(m,n) \f$ singular values, in decreasing order.
     !! @param[in] overwrite_a (Optional) If true, `a` may be overwritten and destroyed. Default is false.
     !! @param[out] rank (Optional) The rank of the matrix `a`.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     public :: solve_lstsq

     !> @brief Compute a weighted least squares solution of \f$ A \cdot x = b \f$.
     !!
     !! This function minimizes \f$ \|D (b - A \cdot x)\| \f$ with \f$ D = \mathrm{diag}(\sqrt{w}) \f$,
     !! that is, it gives the \f$i\f$-th equation the weight \f$ w_i \f$.
     !!
     !! @param[in] w The weight vector of size \f$m\f$. It is always real, and every entry must be positive.
     !! @param[in,out] a The input matrix of size \f$ [m,n] \f$. If `overwrite_a` is true,
     !!                  the contents of `a` may be modified during computation.
     !! @param[in] b The right-hand side vector of size \f$m\f$.
     !! @param[in] cond (Optional) A cutoff for rank evaluation.
     !! @param[in] overwrite_a (Optional) If true, `a` may be overwritten and destroyed. Default is false.
     !! @param[out] rank (Optional) The rank of the matrix `a`.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     !! @return Solution vector `x` of size \f$n\f$.
     !!
     public :: weighted_lstsq

     !> @brief Compute a weighted least squares solution into a pre-allocated array.
     !!
     !! This subroutine is the subroutine form of [`weighted_lstsq`](@ref la_least_squares::weighted_lstsq):
     !! it writes the solution into the caller's array.
     !!
     !! @param[in] w The weight vector of size \f$m\f$. It is always real, and every entry must be positive.
     !! @param[in,out] a The input matrix of size \f$ [m,n] \f$. If `overwrite_a` is true,
     !!                  the contents of `a` may be modified during computation.
     !! @param[in] b The right-hand side vector of size \f$m\f$.
     !! @param[in,out] x The solution vector of size \f$n\f$.
     !! @param[in] cond (Optional) A cutoff for rank evaluation.
     !! @param[in] overwrite_a (Optional) If true, `a` may be overwritten and destroyed. Default is false.
     !! @param[out] rank (Optional) The rank of the matrix `a`.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     public :: solve_weighted_lstsq

     !> @brief Compute the solution of an equality-constrained least squares problem.
     !!
     !! This function minimizes \f$ \|b - A \cdot x\| \f$ subject to \f$ C \cdot x = d \f$,
     !! with \f$ A \f$ of size \f$ [m,n] \f$ and \f$ C \f$ of size \f$ [p,n] \f$,
     !! \f$ p \le n \le m+p \f$.
     !!
     !! @param[in,out] A The least-squares matrix of size \f$ [m,n] \f$.
     !! @param[in,out] b The least-squares right-hand side vector of size \f$m\f$.
     !! @param[in,out] C The constraint matrix of size \f$ [p,n] \f$.
     !! @param[in,out] d The constraint right-hand side vector of size \f$p\f$.
     !! @param[in] overwrite_matrices (Optional) If true, `A`, `b`, `C` and `d` may be overwritten
     !!                               and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     !! @return Solution vector `x` of size \f$n\f$.
     !!
     !! @note This function relies on LAPACK's [GGLSE](@ref la_lapack::gglse) driver.
     !!
     public :: constrained_lstsq

     !> @brief Compute the solution of an equality-constrained least squares problem into a
     !!        pre-allocated array.
     !!
     !! This subroutine is the subroutine form of
     !! [`constrained_lstsq`](@ref la_least_squares::constrained_lstsq): it writes the solution into
     !! the caller's array and can reuse a caller-provided workspace.
     !!
     !! @param[in,out] A The least-squares matrix of size \f$ [m,n] \f$.
     !! @param[in,out] b The least-squares right-hand side vector of size \f$m\f$.
     !! @param[in,out] C The constraint matrix of size \f$ [p,n] \f$.
     !! @param[in,out] d The constraint right-hand side vector of size \f$p\f$.
     !! @param[out] x The solution vector of size \f$n\f$.
     !! @param[out] storage (Optional) Pre-allocated workspace. Its size is checked with
     !!             [`constrained_lstsq_space`](@ref la_least_squares::constrained_lstsq_space).
     !! @param[in] overwrite_matrices (Optional) If true, `A`, `b`, `C` and `d` may be overwritten
     !!                               and destroyed. Default is false.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     public :: solve_constrained_lstsq

     !> @brief Compute the working space needed by the equality-constrained least squares solver.
     !!
     !! This subroutine queries LAPACK for the optimal size of the workspace array that
     !! [`solve_constrained_lstsq`](@ref la_least_squares::solve_constrained_lstsq) needs.
     !!
     !! @param[in] A The least-squares matrix of size \f$ [m,n] \f$. Only its shape is used.
     !! @param[in] C The constraint matrix of size \f$ [p,n] \f$. Only its shape is used.
     !! @param[out] lwork The size of the workspace array.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                    the function will stop execution.
     !!
     public :: constrained_lstsq_space

     interface lstsq
        module procedure la_slstsq_one
        module procedure la_dlstsq_one
        module procedure la_qlstsq_one
        module procedure la_clstsq_one
        module procedure la_zlstsq_one
        module procedure la_wlstsq_one
        module procedure la_slstsq_multiple
        module procedure la_dlstsq_multiple
        module procedure la_qlstsq_multiple
        module procedure la_clstsq_multiple
        module procedure la_zlstsq_multiple
        module procedure la_wlstsq_multiple
     end interface lstsq

     interface lstsq_space
        module procedure la_slstsq_space_one
        module procedure la_dlstsq_space_one
        module procedure la_qlstsq_space_one
        module procedure la_clstsq_space_one
        module procedure la_zlstsq_space_one
        module procedure la_wlstsq_space_one
        module procedure la_slstsq_space_multiple
        module procedure la_dlstsq_space_multiple
        module procedure la_qlstsq_space_multiple
        module procedure la_clstsq_space_multiple
        module procedure la_zlstsq_space_multiple
        module procedure la_wlstsq_space_multiple
     end interface lstsq_space

     interface solve_lstsq
        module procedure la_ssolve_lstsq_one
        module procedure la_dsolve_lstsq_one
        module procedure la_qsolve_lstsq_one
        module procedure la_csolve_lstsq_one
        module procedure la_zsolve_lstsq_one
        module procedure la_wsolve_lstsq_one
        module procedure la_ssolve_lstsq_multiple
        module procedure la_dsolve_lstsq_multiple
        module procedure la_qsolve_lstsq_multiple
        module procedure la_csolve_lstsq_multiple
        module procedure la_zsolve_lstsq_multiple
        module procedure la_wsolve_lstsq_multiple
     end interface solve_lstsq

     interface weighted_lstsq
        module procedure la_sweighted_lstsq
        module procedure la_dweighted_lstsq
        module procedure la_qweighted_lstsq
        module procedure la_cweighted_lstsq
        module procedure la_zweighted_lstsq
        module procedure la_wweighted_lstsq
     end interface weighted_lstsq

     interface solve_weighted_lstsq
        module procedure la_ssolve_weighted_lstsq
        module procedure la_dsolve_weighted_lstsq
        module procedure la_qsolve_weighted_lstsq
        module procedure la_csolve_weighted_lstsq
        module procedure la_zsolve_weighted_lstsq
        module procedure la_wsolve_weighted_lstsq
     end interface solve_weighted_lstsq

     interface constrained_lstsq
        module procedure la_sconstrained_lstsq
        module procedure la_dconstrained_lstsq
        module procedure la_qconstrained_lstsq
        module procedure la_cconstrained_lstsq
        module procedure la_zconstrained_lstsq
        module procedure la_wconstrained_lstsq
     end interface constrained_lstsq

     interface solve_constrained_lstsq
        module procedure la_ssolve_constrained_lstsq
        module procedure la_dsolve_constrained_lstsq
        module procedure la_qsolve_constrained_lstsq
        module procedure la_csolve_constrained_lstsq
        module procedure la_zsolve_constrained_lstsq
        module procedure la_wsolve_constrained_lstsq
     end interface solve_constrained_lstsq

     interface constrained_lstsq_space
        module procedure la_sconstrained_lstsq_space
        module procedure la_dconstrained_lstsq_space
        module procedure la_qconstrained_lstsq_space
        module procedure la_cconstrained_lstsq_space
        module procedure la_zconstrained_lstsq_space
        module procedure la_wconstrained_lstsq_space
     end interface constrained_lstsq_space

     contains

     !> Workspace needed by real(sp) gelsd
     pure subroutine sgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'sgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 12*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + mnmin*nrhs + (smlsiz + 1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine sgelsd_space

     !> Workspace needed by real(dp) gelsd
     pure subroutine dgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'dgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 12*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + mnmin*nrhs + (smlsiz + 1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine dgelsd_space

     !> Workspace needed by real(qp) gelsd
     pure subroutine qgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'qgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 12*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + mnmin*nrhs + (smlsiz + 1)**2
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine qgelsd_space

     !> Workspace needed by complex(sp) gelsd
     pure subroutine cgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'cgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 10*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*(1 + nrhs) + 2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine cgelsd_space

     !> Workspace needed by complex(dp) gelsd
     pure subroutine zgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'zgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 10*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*(1 + nrhs) + 2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine zgelsd_space

     !> Workspace needed by complex(qp) gelsd
     pure subroutine wgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)
         integer(ilp),intent(in) :: m,n,nrhs
         integer(ilp),intent(out) :: lrwork,liwork,lcwork

         integer(ilp) :: smlsiz,mnmin,nlvl

         mnmin = min(m,n)

         ! Maximum size of the subproblems at the bottom of the computation (~25)
         smlsiz = la_ilaenv(9,'wgelsd',' ',0,0,0,0)

         ! The exact minimum amount of workspace needed depends on M, N and NRHS. As long as LWORK is at least
         nlvl = max(0,ilog2(mnmin/(smlsiz + 1)) + 1)

         ! Real space
         lrwork = 10*mnmin + 2*mnmin*smlsiz + 8*mnmin*nlvl + 3*smlsiz*nrhs + max((smlsiz + 1)**2,n*(1 + nrhs) + 2*nrhs)
         lrwork = max(1,lrwork)

         ! Complex space
         lcwork = 2*mnmin + nrhs*mnmin

         ! Integer space
         liwork = max(1,3*mnmin*nlvl + 11*mnmin)

         ! For good performance, the workspace should generally be larger.
         lrwork = ceiling(1.25*lrwork,kind=ilp)
         lcwork = ceiling(1.25*lcwork,kind=ilp)
         liwork = ceiling(1.25*liwork,kind=ilp)

     end subroutine wgelsd_space

     !> Working space needed by the real(sp) least-squares solver
     pure subroutine la_slstsq_space_one(a,b,lrwork,liwork)
         !> Input matrix a[m,n]
         real(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         real(sp),intent(in) :: b(:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         integer(ilp) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call sgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_slstsq_space_one

     !> Compute the least-squares solution to a real(sp) system of linear equations Ax = B
     function la_slstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n))

         call la_ssolve_lstsq_one(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0))
         end if

         call err0%handle(err)

     end function la_slstsq_one

     !> Compute the least-squares solution to a real(sp) system of linear equations Ax = B into x
     subroutine la_ssolve_lstsq_one(a,b,x,real_storage,int_storage, &
                cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] real working storage space
         real(sp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(sp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(sp) :: rcond
         real(sp),pointer :: rwork(:),singular(:)
         real(sp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Get working space size
         call sgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       rwork,nrs,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n) = xmat(1:n,1)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_lstsq_one

     !> Working space needed by the real(dp) least-squares solver
     pure subroutine la_dlstsq_space_one(a,b,lrwork,liwork)
         !> Input matrix a[m,n]
         real(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         real(dp),intent(in) :: b(:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         integer(ilp) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call dgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_dlstsq_space_one

     !> Compute the least-squares solution to a real(dp) system of linear equations Ax = B
     function la_dlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n))

         call la_dsolve_lstsq_one(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0))
         end if

         call err0%handle(err)

     end function la_dlstsq_one

     !> Compute the least-squares solution to a real(dp) system of linear equations Ax = B into x
     subroutine la_dsolve_lstsq_one(a,b,x,real_storage,int_storage, &
                cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] real working storage space
         real(dp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(dp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(dp) :: rcond
         real(dp),pointer :: rwork(:),singular(:)
         real(dp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Get working space size
         call dgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       rwork,nrs,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n) = xmat(1:n,1)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_lstsq_one

     !> Working space needed by the real(qp) least-squares solver
     pure subroutine la_qlstsq_space_one(a,b,lrwork,liwork)
         !> Input matrix a[m,n]
         real(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         real(qp),intent(in) :: b(:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         integer(ilp) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call qgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_qlstsq_space_one

     !> Compute the least-squares solution to a real(qp) system of linear equations Ax = B
     function la_qlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n))

         call la_qsolve_lstsq_one(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0))
         end if

         call err0%handle(err)

     end function la_qlstsq_one

     !> Compute the least-squares solution to a real(qp) system of linear equations Ax = B into x
     subroutine la_qsolve_lstsq_one(a,b,x,real_storage,int_storage, &
                cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] real working storage space
         real(qp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(qp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(qp) :: rcond
         real(qp),pointer :: rwork(:),singular(:)
         real(qp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Get working space size
         call qgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       rwork,nrs,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n) = xmat(1:n,1)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_lstsq_one

     !> Working space needed by the complex(sp) least-squares solver
     pure subroutine la_clstsq_space_one(a,b,lrwork,liwork,lcwork)
         !> Input matrix a[m,n]
         complex(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         !> Size of the complex working space array
         integer(ilp),intent(out) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call cgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_clstsq_space_one

     !> Compute the least-squares solution to a complex(sp) system of linear equations Ax = B
     function la_clstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n))

         call la_csolve_lstsq_one(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0))
         end if

         call err0%handle(err)

     end function la_clstsq_one

     !> Compute the least-squares solution to a complex(sp) system of linear equations Ax = B into x
     subroutine la_csolve_lstsq_one(a,b,x,real_storage,int_storage, &
                cmpl_storage,cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] real working storage space
         real(sp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] complex working storage space
         complex(sp),optional,intent(inout),target :: cmpl_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(sp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp) :: ncs
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(sp) :: rcond
         real(sp),pointer :: rwork(:),singular(:)
         complex(sp),pointer :: xmat(:,:),amat(:,:)
         complex(sp),pointer :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Get working space size
         call cgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         ! Complex working space
         if (present(cmpl_storage)) then
            cwork => cmpl_storage
         else
            allocate (cwork(lcwork))
         end if
         ncs = size(cwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork .or. ncs < lcwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', cmplx=',ncs,' should be >=',lcwork, &
                       ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       cwork,ncs,rwork,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n) = xmat(1:n,1)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(cmpl_storage)) deallocate (cwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_lstsq_one

     !> Working space needed by the complex(dp) least-squares solver
     pure subroutine la_zlstsq_space_one(a,b,lrwork,liwork,lcwork)
         !> Input matrix a[m,n]
         complex(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         !> Size of the complex working space array
         integer(ilp),intent(out) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call zgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_zlstsq_space_one

     !> Compute the least-squares solution to a complex(dp) system of linear equations Ax = B
     function la_zlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n))

         call la_zsolve_lstsq_one(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0))
         end if

         call err0%handle(err)

     end function la_zlstsq_one

     !> Compute the least-squares solution to a complex(dp) system of linear equations Ax = B into x
     subroutine la_zsolve_lstsq_one(a,b,x,real_storage,int_storage, &
                cmpl_storage,cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] real working storage space
         real(dp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] complex working storage space
         complex(dp),optional,intent(inout),target :: cmpl_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(dp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp) :: ncs
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(dp) :: rcond
         real(dp),pointer :: rwork(:),singular(:)
         complex(dp),pointer :: xmat(:,:),amat(:,:)
         complex(dp),pointer :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Get working space size
         call zgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         ! Complex working space
         if (present(cmpl_storage)) then
            cwork => cmpl_storage
         else
            allocate (cwork(lcwork))
         end if
         ncs = size(cwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork .or. ncs < lcwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', cmplx=',ncs,' should be >=',lcwork, &
                       ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       cwork,ncs,rwork,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n) = xmat(1:n,1)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(cmpl_storage)) deallocate (cwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_lstsq_one

     !> Working space needed by the complex(qp) least-squares solver
     pure subroutine la_wlstsq_space_one(a,b,lrwork,liwork,lcwork)
         !> Input matrix a[m,n]
         complex(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         !> Size of the complex working space array
         integer(ilp),intent(out) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call wgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_wlstsq_space_one

     !> Compute the least-squares solution to a complex(qp) system of linear equations Ax = B
     function la_wlstsq_one(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n))

         call la_wsolve_lstsq_one(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0))
         end if

         call err0%handle(err)

     end function la_wlstsq_one

     !> Compute the least-squares solution to a complex(qp) system of linear equations Ax = B into x
     subroutine la_wsolve_lstsq_one(a,b,x,real_storage,int_storage, &
                cmpl_storage,cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] real working storage space
         real(qp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] complex working storage space
         complex(qp),optional,intent(inout),target :: cmpl_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(qp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp) :: ncs
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(qp) :: rcond
         real(qp),pointer :: rwork(:),singular(:)
         complex(qp),pointer :: xmat(:,:),amat(:,:)
         complex(qp),pointer :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Get working space size
         call wgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         ! Complex working space
         if (present(cmpl_storage)) then
            cwork => cmpl_storage
         else
            allocate (cwork(lcwork))
         end if
         ncs = size(cwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork .or. ncs < lcwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', cmplx=',ncs,' should be >=',lcwork, &
                       ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       cwork,ncs,rwork,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n) = xmat(1:n,1)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(cmpl_storage)) deallocate (cwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_lstsq_one

     !> Working space needed by the real(sp) least-squares solver
     pure subroutine la_slstsq_space_multiple(a,b,lrwork,liwork)
         !> Input matrix a[m,n]
         real(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         integer(ilp) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call sgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_slstsq_space_multiple

     !> Compute the least-squares solution to a real(sp) system of linear equations Ax = B
     function la_slstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n,nrhs))

         call la_ssolve_lstsq_multiple(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0,0))
         end if

         call err0%handle(err)

     end function la_slstsq_multiple

     !> Compute the least-squares solution to a real(sp) system of linear equations Ax = B into x
     subroutine la_ssolve_lstsq_multiple(a,b,x,real_storage,int_storage, &
                cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] real working storage space
         real(sp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(sp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(sp) :: rcond
         real(sp),pointer :: rwork(:),singular(:)
         real(sp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1:nrhs) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Get working space size
         call sgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       rwork,nrs,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n,1:nrhs) = xmat(1:n,1:nrhs)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_lstsq_multiple

     !> Working space needed by the real(dp) least-squares solver
     pure subroutine la_dlstsq_space_multiple(a,b,lrwork,liwork)
         !> Input matrix a[m,n]
         real(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         integer(ilp) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call dgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_dlstsq_space_multiple

     !> Compute the least-squares solution to a real(dp) system of linear equations Ax = B
     function la_dlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n,nrhs))

         call la_dsolve_lstsq_multiple(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0,0))
         end if

         call err0%handle(err)

     end function la_dlstsq_multiple

     !> Compute the least-squares solution to a real(dp) system of linear equations Ax = B into x
     subroutine la_dsolve_lstsq_multiple(a,b,x,real_storage,int_storage, &
                cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] real working storage space
         real(dp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(dp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(dp) :: rcond
         real(dp),pointer :: rwork(:),singular(:)
         real(dp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1:nrhs) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Get working space size
         call dgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       rwork,nrs,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n,1:nrhs) = xmat(1:n,1:nrhs)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_lstsq_multiple

     !> Working space needed by the real(qp) least-squares solver
     pure subroutine la_qlstsq_space_multiple(a,b,lrwork,liwork)
         !> Input matrix a[m,n]
         real(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         integer(ilp) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call qgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_qlstsq_space_multiple

     !> Compute the least-squares solution to a real(qp) system of linear equations Ax = B
     function la_qlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n,nrhs))

         call la_qsolve_lstsq_multiple(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0,0))
         end if

         call err0%handle(err)

     end function la_qlstsq_multiple

     !> Compute the least-squares solution to a real(qp) system of linear equations Ax = B into x
     subroutine la_qsolve_lstsq_multiple(a,b,x,real_storage,int_storage, &
                cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         real(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         real(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] real working storage space
         real(qp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(qp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(qp) :: rcond
         real(qp),pointer :: rwork(:),singular(:)
         real(qp),pointer :: xmat(:,:),amat(:,:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1:nrhs) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Get working space size
         call qgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       rwork,nrs,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n,1:nrhs) = xmat(1:n,1:nrhs)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_lstsq_multiple

     !> Working space needed by the complex(sp) least-squares solver
     pure subroutine la_clstsq_space_multiple(a,b,lrwork,liwork,lcwork)
         !> Input matrix a[m,n]
         complex(sp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         !> Size of the complex working space array
         integer(ilp),intent(out) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call cgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_clstsq_space_multiple

     !> Compute the least-squares solution to a complex(sp) system of linear equations Ax = B
     function la_clstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n,nrhs))

         call la_csolve_lstsq_multiple(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0,0))
         end if

         call err0%handle(err)

     end function la_clstsq_multiple

     !> Compute the least-squares solution to a complex(sp) system of linear equations Ax = B into x
     subroutine la_csolve_lstsq_multiple(a,b,x,real_storage,int_storage, &
                cmpl_storage,cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(sp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(sp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] real working storage space
         real(sp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] complex working storage space
         complex(sp),optional,intent(inout),target :: cmpl_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(sp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp) :: ncs
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(sp) :: rcond
         real(sp),pointer :: rwork(:),singular(:)
         complex(sp),pointer :: xmat(:,:),amat(:,:)
         complex(sp),pointer :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1:nrhs) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_sp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_sp)*mnmax

         ! Get working space size
         call cgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         ! Complex working space
         if (present(cmpl_storage)) then
            cwork => cmpl_storage
         else
            allocate (cwork(lcwork))
         end if
         ncs = size(cwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork .or. ncs < lcwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', cmplx=',ncs,' should be >=',lcwork, &
                       ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       cwork,ncs,rwork,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n,1:nrhs) = xmat(1:n,1:nrhs)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(cmpl_storage)) deallocate (cwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_lstsq_multiple

     !> Working space needed by the complex(dp) least-squares solver
     pure subroutine la_zlstsq_space_multiple(a,b,lrwork,liwork,lcwork)
         !> Input matrix a[m,n]
         complex(dp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         !> Size of the complex working space array
         integer(ilp),intent(out) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call zgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_zlstsq_space_multiple

     !> Compute the least-squares solution to a complex(dp) system of linear equations Ax = B
     function la_zlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n,nrhs))

         call la_zsolve_lstsq_multiple(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0,0))
         end if

         call err0%handle(err)

     end function la_zlstsq_multiple

     !> Compute the least-squares solution to a complex(dp) system of linear equations Ax = B into x
     subroutine la_zsolve_lstsq_multiple(a,b,x,real_storage,int_storage, &
                cmpl_storage,cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(dp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(dp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] real working storage space
         real(dp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] complex working storage space
         complex(dp),optional,intent(inout),target :: cmpl_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(dp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp) :: ncs
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(dp) :: rcond
         real(dp),pointer :: rwork(:),singular(:)
         complex(dp),pointer :: xmat(:,:),amat(:,:)
         complex(dp),pointer :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1:nrhs) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_dp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_dp)*mnmax

         ! Get working space size
         call zgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         ! Complex working space
         if (present(cmpl_storage)) then
            cwork => cmpl_storage
         else
            allocate (cwork(lcwork))
         end if
         ncs = size(cwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork .or. ncs < lcwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', cmplx=',ncs,' should be >=',lcwork, &
                       ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       cwork,ncs,rwork,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n,1:nrhs) = xmat(1:n,1:nrhs)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(cmpl_storage)) deallocate (cwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_lstsq_multiple

     !> Working space needed by the complex(qp) least-squares solver
     pure subroutine la_wlstsq_space_multiple(a,b,lrwork,liwork,lcwork)
         !> Input matrix a[m,n]
         complex(qp),intent(in) :: a(:,:)
         !> Right hand side vector or array, b[m] or b[m,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Size of the real working space array
         integer(ilp),intent(out) :: lrwork
         !> Size of the integer working space array
         integer(ilp),intent(out) :: liwork
         !> Size of the complex working space array
         integer(ilp),intent(out) :: lcwork

         integer(ilp) :: m,n,nrhs

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         nrhs = size(b,kind=ilp)/size(b,1,kind=ilp)

         call wgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

     end subroutine la_wlstsq_space_multiple

     !> Compute the least-squares solution to a complex(qp) system of linear equations Ax = B
     function la_wlstsq_multiple(a,b,cond,overwrite_a,rank,err) result(x)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),allocatable,target :: x(:,:)

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: n,nrhs,ldb

         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb

         ! The solution has the shape of the second dimension of a
         allocate (x(n,nrhs))

         call la_wsolve_lstsq_multiple(a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err0)

         ! An invalid problem returns an empty solution
         if (err0%error()) then
            deallocate (x)
            allocate (x(0,0))
         end if

         call err0%handle(err)

     end function la_wlstsq_multiple

     !> Compute the least-squares solution to a complex(qp) system of linear equations Ax = B into x
     subroutine la_wsolve_lstsq_multiple(a,b,x,real_storage,int_storage, &
                cmpl_storage,cond,singvals,overwrite_a,rank,err)
         !> Input matrix a[n,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector or array, b[n] or b[n,nrhs]
         complex(qp),intent(in) :: b(:,:)
         !> Result array/matrix x[n] or x[n,nrhs]
         complex(qp),intent(inout),contiguous,target :: x(:,:)
         !> [optional] real working storage space
         real(qp),optional,intent(inout),target :: real_storage(:)
         !> [optional] integer working storage space
         integer(ilp),optional,intent(inout),target :: int_storage(:)
         !> [optional] complex working storage space
         complex(qp),optional,intent(inout),target :: cmpl_storage(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] list of singular values [min(m,n)], in descending magnitude order, returned by the SVD
         real(qp),optional,intent(out),target :: singvals(:)
         !> [optional] Can A,b data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,lda,ldb,ldx,nrhs,nrhsx,info,mnmin,mnmax,arank,lrwork,liwork,lcwork
         integer(ilp) :: nrs,nis,nsvd
         integer(ilp) :: ncs
         integer(ilp),pointer :: iwork(:)
         logical(lk) :: copy_a,large_enough_x
         real(qp) :: rcond
         real(qp),pointer :: rwork(:),singular(:)
         complex(qp),pointer :: xmat(:,:),amat(:,:)
         complex(qp),pointer :: cwork(:)
         character(*),parameter :: this = 'lstsq'

         !> Problem sizes
         m = size(a,1,kind=ilp)
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         ldb = size(b,1,kind=ilp)
         nrhs = size(b,kind=ilp)/ldb
         ldx = size(x,1,kind=ilp)
         nrhsx = size(x,kind=ilp)/ldx
         mnmin = min(m,n)
         mnmax = max(m,n)
         arank = 0

         if (lda < 1 .or. n < 1 .or. ldb < 1 .or. ldb /= m .or. ldx < n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient sizes: a=[',lda,',',n,'],', &
                                                                    'b=[',ldb,',',nrhs,'],', &
                                                                    'x=[',ldx,',',nrhsx,']')
            call err0%handle(err)
            if (present(rank)) rank = arank
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

         ! *GELSD stores the rhs, and returns the solution, in an array with max(m,n) rows:
         ! x can host it whenever its leading dimension is at least m, otherwise a temporary is needed
         large_enough_x = ldx >= m
         if (large_enough_x) then
            xmat(1:ldx,1:nrhs) => x
         else
            allocate (xmat(m,nrhs))
         end if
         xmat(1:m,1:nrhs) = b

         ! Singular values array (in decreasing order)
         if (present(singvals)) then
            singular => singvals
            nsvd = size(singular,kind=ilp)
         else
            allocate (singular(mnmin))
            nsvd = mnmin
         end if

         ! rcond is used to determine the effective rank of A.
         ! Singular values S(i) <= RCOND*maxval(S) are treated as zero.
         ! Use same default value as NumPy
         if (present(cond)) then
            rcond = cond
         else
            rcond = epsilon(0.0_qp)*mnmax
         end if
         if (rcond < 0) rcond = epsilon(0.0_qp)*mnmax

         ! Get working space size
         call wgelsd_space(m,n,nrhs,lrwork,liwork,lcwork)

         ! Real working space
         if (present(real_storage)) then
            rwork => real_storage
         else
            allocate (rwork(lrwork))
         end if
         nrs = size(rwork,kind=ilp)

         ! Integer working space
         if (present(int_storage)) then
            iwork => int_storage
         else
            allocate (iwork(liwork))
         end if
         nis = size(iwork,kind=ilp)

         ! Complex working space
         if (present(cmpl_storage)) then
            cwork => cmpl_storage
         else
            allocate (cwork(lcwork))
         end if
         ncs = size(cwork,kind=ilp)

         if (nrs < lrwork .or. nis < liwork .or. ncs < lcwork &
             .or. nsvd < mnmin) then

            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:', &
                          'real=',nrs,' should be >=',lrwork, &
                         ', int=',nis,' should be >=',liwork, &
                         ', cmplx=',ncs,' should be >=',lcwork, &
                       ', singv=',nsvd,' should be >=',mnmin)

         else

            ! Solve system using singular value decomposition
            call gelsd(m,n,nrhs,amat,lda,xmat,size(xmat,1,kind=ilp),singular,rcond,arank, &
                       cwork,ncs,rwork,iwork,info)

            ! Process output
            select case (info)
               case (0)
                   ! Success
               case (:-1)
                   err0 = la_state(this,LINALG_VALUE_ERROR,'invalid problem size a=', [lda,n],', b=', [ldb,nrhs])
               case (1:)
                   err0 = la_state(this,LINALG_ERROR,'SVD did not converge.')
               case default
                   err0 = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select

         end if

         ! Retrieve the solution from the temporary storage
         if (.not. large_enough_x) then
            x(1:n,1:nrhs) = xmat(1:n,1:nrhs)
            deallocate (xmat)
         end if

         if (copy_a) deallocate (amat)
         if (present(rank)) rank = arank
         if (.not. present(real_storage)) deallocate (rwork)
         if (.not. present(int_storage)) deallocate (iwork)
         if (.not. present(cmpl_storage)) deallocate (cwork)
         if (.not. present(singvals)) deallocate (singular)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_lstsq_multiple

     ! Simple integer log2 implementation
     elemental integer(ilp) function ilog2(x)
        integer(ilp),intent(in) :: x

        integer(ilp) :: remndr

        if (x > 0) then
           remndr = x
           ilog2 = -1_ilp
           do while (remndr > 0)
               ilog2 = ilog2 + 1_ilp
               remndr = shiftr(remndr,1)
           end do
        else
           ilog2 = -huge(0_ilp)
        end if
     end function ilog2

     !> Validate the sizes of an equality-constrained least-squares problem
     pure subroutine la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err)
         integer(ilp),intent(in) :: ma,na,mb,mc,nc,md,mx
         type(la_state),intent(out) :: err

         character(*),parameter :: this = 'constrained_lstsq'

         if (ma < 1 .or. na < 1) then
            err = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [ma,na])
         elseif (mc < 1 .or. nc < 1) then
            err = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size c(p,n) =', [mc,nc])
         elseif (na /= nc) then
            err = la_state(this,LINALG_VALUE_ERROR,'a and c have',na,'and',nc,'columns')
         elseif (mb /= ma) then
            err = la_state(this,LINALG_VALUE_ERROR,'size(b) =',mb,'but a has',ma,'rows')
         elseif (md /= mc) then
            err = la_state(this,LINALG_VALUE_ERROR,'size(d) =',md,'but c has',mc,'rows')
         elseif (na /= mx) then
            err = la_state(this,LINALG_VALUE_ERROR,'size(x) =',mx,'but a has',na,'columns')
         end if

     end subroutine la_check_constrained_sizes

     !> Process the gglse exit code
     pure subroutine la_handle_gglse(info,m,n,p,err)
         integer(ilp),intent(in) :: info,m,n,p
         type(la_state),intent(out) :: err

         character(*),parameter :: this = 'constrained_lstsq'

         select case (info)
            case (0)
                ! Success
            case (1)
                err = la_state(this,LINALG_ERROR,'rank(c) < p: the least-squares solution cannot be computed')
            case (2)
                err = la_state(this,LINALG_ERROR,'rank([a; c]) < n: the least-squares solution cannot be computed')
            case (-1)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid number of rows for a, m=',m)
            case (-2)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid number of columns for a and c, n=',n)
            case (-3)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid number of rows for c, p=',p)
            case default
                err = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine la_handle_gglse

     !> Working space needed by the real(sp) equality-constrained least-squares solver
     subroutine la_sconstrained_lstsq_space(A,C,lwork,err)
         !> Least-squares matrix a[m,n]
         real(sp),intent(in) :: A(:,:)
         !> Equality constraint matrix c[p,n]
         real(sp),intent(in) :: C(:,:)
         !> Size of the working space array
         integer(ilp),intent(out) :: lwork
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,p,info
         real(sp) :: a_dummy(1,1),b_dummy(1)
         real(sp) :: c_dummy(1,1),d_dummy(1)
         real(sp) :: work(1),x(1)

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         p = size(C,1,kind=ilp)

         lwork = -1_ilp
         call gglse(m,n,p,a_dummy,m,c_dummy,p,b_dummy,d_dummy,x,work,lwork,info)
         call la_handle_gglse(info,m,n,p,err0)

         lwork = ceiling(real(work(1),kind=sp),kind=ilp)

         call err0%handle(err)

     end subroutine la_sconstrained_lstsq_space

     !> Solve the real(sp) equality-constrained least-squares problem into x
     subroutine la_ssolve_constrained_lstsq(A,b,C,d,x,storage,overwrite_matrices,err)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         real(sp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         real(sp),intent(inout),target :: C(:,:),d(:)
         !> Solution vector x[n]
         real(sp),intent(out) :: x(:)
         !> [optional] working storage space
         real(sp),optional,intent(out),target :: storage(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: ma,na,mb,mc,nc,md,mx,lwork,info
         logical(lk) :: copy_matrices
         real(sp),pointer :: amat(:,:),bvec(:),cmat(:,:),dvec(:),work(:)
         character(*),parameter :: this = 'constrained_lstsq'

         ma = size(A,1,kind=ilp)
         na = size(A,2,kind=ilp)
         mc = size(C,1,kind=ilp)
         nc = size(C,2,kind=ilp)
         mb = size(b,kind=ilp)
         md = size(d,kind=ilp)
         mx = size(x,kind=ilp)

         call la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err0)
         if (err0%error()) then
            call err0%handle(err)
            return
         end if

         ! Can the input matrices be overwritten? By default, do not overwrite
         if (present(overwrite_matrices)) then
            copy_matrices = .not. overwrite_matrices
         else
            copy_matrices = .true._lk
         end if

         if (copy_matrices) then
            allocate (amat(ma,na),source=A)
            allocate (bvec(mb),source=b)
            allocate (cmat(mc,nc),source=C)
            allocate (dvec(md),source=d)
         else
            amat => A
            bvec => b
            cmat => C
            dvec => d
         end if

         call la_sconstrained_lstsq_space(A,C,lwork,err0)

         if (err0%ok()) then

            if (present(storage)) then
               work => storage
            else
               allocate (work(lwork))
            end if

            if (size(work,kind=ilp) < lwork) then
               err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:',size(work,kind=ilp), &
                                                       ' should be >=',lwork)
            else
               call gglse(ma,na,mc,amat,ma,cmat,mc,bvec,dvec,x,work,lwork,info)
               call la_handle_gglse(info,ma,na,mc,err0)
            end if

            if (.not. present(storage)) deallocate (work)

         end if

         if (copy_matrices) deallocate (amat,bvec,cmat,dvec)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_constrained_lstsq

     !> Solve the real(sp) equality-constrained least-squares problem
     function la_sconstrained_lstsq(A,b,C,d,overwrite_matrices,err) result(x)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         real(sp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         real(sp),intent(inout),target :: C(:,:),d(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Solution vector x[n]
         real(sp),allocatable,target :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(A,2,kind=ilp)
         allocate (x(n))

         call la_ssolve_constrained_lstsq(A,b,C,d,x,overwrite_matrices=overwrite_matrices,err=err)

     end function la_sconstrained_lstsq

     !> Compute the real(sp) weighted least-squares solution
     function la_sweighted_lstsq(w,a,b,cond,overwrite_a,rank,err) result(x)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(sp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         real(sp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array x[n]
         real(sp),allocatable :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(a,2,kind=ilp)
         allocate (x(n))

         call la_ssolve_weighted_lstsq(w,a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err)

     end function la_sweighted_lstsq

     !> Compute the real(sp) weighted least-squares solution into x: minimize ||D(Ax - b)||, D = diag(sqrt(w))
     subroutine la_ssolve_weighted_lstsq(w,a,b,x,cond,overwrite_a,rank,err)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(sp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         real(sp),intent(in) :: b(:)
         !> Result array x[n]
         real(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,i
         logical(lk) :: copy_a
         real(sp),pointer :: amat(:,:)
         real(sp),allocatable,target :: amat_alloc(:,:)
         real(sp),allocatable :: b_scaled(:)
         real(sp),allocatable :: sqrt_w(:)
         character(*),parameter :: this = 'weighted_lstsq'

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m < 1 .or. n < 1) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [m,n])
         elseif (size(w,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(w) =',size(w,kind=ilp),'but a has',m,'rows')
         elseif (size(b,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(b) =',size(b,kind=ilp),'but a has',m,'rows')
         elseif (any(w <= 0.0_sp)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'the weights must be positive')
         end if

         if (err0%error()) then
            call err0%handle(err)
            if (present(rank)) rank = 0_ilp
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (copy_a) then
            allocate (amat_alloc(m,n),source=a)
            amat => amat_alloc
         else
            amat => a
         end if

         ! Scale the rows of A and of b by the square root of the weights
         sqrt_w = sqrt(w)
         do i = 1,m
            amat(i,:) = sqrt_w(i)*amat(i,:)
         end do
         b_scaled = sqrt_w*b

         call la_ssolve_lstsq_one(amat,b_scaled,x,cond=cond,overwrite_a=.true._lk,rank=rank,err=err0)

         ! Report the error against this procedure, not against the transformed problem
         if (err0%error()) err0%where_at = this

         if (copy_a) deallocate (amat_alloc)
         deallocate (b_scaled,sqrt_w)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_ssolve_weighted_lstsq

     !> Working space needed by the real(dp) equality-constrained least-squares solver
     subroutine la_dconstrained_lstsq_space(A,C,lwork,err)
         !> Least-squares matrix a[m,n]
         real(dp),intent(in) :: A(:,:)
         !> Equality constraint matrix c[p,n]
         real(dp),intent(in) :: C(:,:)
         !> Size of the working space array
         integer(ilp),intent(out) :: lwork
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,p,info
         real(dp) :: a_dummy(1,1),b_dummy(1)
         real(dp) :: c_dummy(1,1),d_dummy(1)
         real(dp) :: work(1),x(1)

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         p = size(C,1,kind=ilp)

         lwork = -1_ilp
         call gglse(m,n,p,a_dummy,m,c_dummy,p,b_dummy,d_dummy,x,work,lwork,info)
         call la_handle_gglse(info,m,n,p,err0)

         lwork = ceiling(real(work(1),kind=dp),kind=ilp)

         call err0%handle(err)

     end subroutine la_dconstrained_lstsq_space

     !> Solve the real(dp) equality-constrained least-squares problem into x
     subroutine la_dsolve_constrained_lstsq(A,b,C,d,x,storage,overwrite_matrices,err)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         real(dp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         real(dp),intent(inout),target :: C(:,:),d(:)
         !> Solution vector x[n]
         real(dp),intent(out) :: x(:)
         !> [optional] working storage space
         real(dp),optional,intent(out),target :: storage(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: ma,na,mb,mc,nc,md,mx,lwork,info
         logical(lk) :: copy_matrices
         real(dp),pointer :: amat(:,:),bvec(:),cmat(:,:),dvec(:),work(:)
         character(*),parameter :: this = 'constrained_lstsq'

         ma = size(A,1,kind=ilp)
         na = size(A,2,kind=ilp)
         mc = size(C,1,kind=ilp)
         nc = size(C,2,kind=ilp)
         mb = size(b,kind=ilp)
         md = size(d,kind=ilp)
         mx = size(x,kind=ilp)

         call la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err0)
         if (err0%error()) then
            call err0%handle(err)
            return
         end if

         ! Can the input matrices be overwritten? By default, do not overwrite
         if (present(overwrite_matrices)) then
            copy_matrices = .not. overwrite_matrices
         else
            copy_matrices = .true._lk
         end if

         if (copy_matrices) then
            allocate (amat(ma,na),source=A)
            allocate (bvec(mb),source=b)
            allocate (cmat(mc,nc),source=C)
            allocate (dvec(md),source=d)
         else
            amat => A
            bvec => b
            cmat => C
            dvec => d
         end if

         call la_dconstrained_lstsq_space(A,C,lwork,err0)

         if (err0%ok()) then

            if (present(storage)) then
               work => storage
            else
               allocate (work(lwork))
            end if

            if (size(work,kind=ilp) < lwork) then
               err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:',size(work,kind=ilp), &
                                                       ' should be >=',lwork)
            else
               call gglse(ma,na,mc,amat,ma,cmat,mc,bvec,dvec,x,work,lwork,info)
               call la_handle_gglse(info,ma,na,mc,err0)
            end if

            if (.not. present(storage)) deallocate (work)

         end if

         if (copy_matrices) deallocate (amat,bvec,cmat,dvec)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_constrained_lstsq

     !> Solve the real(dp) equality-constrained least-squares problem
     function la_dconstrained_lstsq(A,b,C,d,overwrite_matrices,err) result(x)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         real(dp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         real(dp),intent(inout),target :: C(:,:),d(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Solution vector x[n]
         real(dp),allocatable,target :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(A,2,kind=ilp)
         allocate (x(n))

         call la_dsolve_constrained_lstsq(A,b,C,d,x,overwrite_matrices=overwrite_matrices,err=err)

     end function la_dconstrained_lstsq

     !> Compute the real(dp) weighted least-squares solution
     function la_dweighted_lstsq(w,a,b,cond,overwrite_a,rank,err) result(x)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(dp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         real(dp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array x[n]
         real(dp),allocatable :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(a,2,kind=ilp)
         allocate (x(n))

         call la_dsolve_weighted_lstsq(w,a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err)

     end function la_dweighted_lstsq

     !> Compute the real(dp) weighted least-squares solution into x: minimize ||D(Ax - b)||, D = diag(sqrt(w))
     subroutine la_dsolve_weighted_lstsq(w,a,b,x,cond,overwrite_a,rank,err)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(dp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         real(dp),intent(in) :: b(:)
         !> Result array x[n]
         real(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,i
         logical(lk) :: copy_a
         real(dp),pointer :: amat(:,:)
         real(dp),allocatable,target :: amat_alloc(:,:)
         real(dp),allocatable :: b_scaled(:)
         real(dp),allocatable :: sqrt_w(:)
         character(*),parameter :: this = 'weighted_lstsq'

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m < 1 .or. n < 1) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [m,n])
         elseif (size(w,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(w) =',size(w,kind=ilp),'but a has',m,'rows')
         elseif (size(b,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(b) =',size(b,kind=ilp),'but a has',m,'rows')
         elseif (any(w <= 0.0_dp)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'the weights must be positive')
         end if

         if (err0%error()) then
            call err0%handle(err)
            if (present(rank)) rank = 0_ilp
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (copy_a) then
            allocate (amat_alloc(m,n),source=a)
            amat => amat_alloc
         else
            amat => a
         end if

         ! Scale the rows of A and of b by the square root of the weights
         sqrt_w = sqrt(w)
         do i = 1,m
            amat(i,:) = sqrt_w(i)*amat(i,:)
         end do
         b_scaled = sqrt_w*b

         call la_dsolve_lstsq_one(amat,b_scaled,x,cond=cond,overwrite_a=.true._lk,rank=rank,err=err0)

         ! Report the error against this procedure, not against the transformed problem
         if (err0%error()) err0%where_at = this

         if (copy_a) deallocate (amat_alloc)
         deallocate (b_scaled,sqrt_w)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_dsolve_weighted_lstsq

     !> Working space needed by the real(qp) equality-constrained least-squares solver
     subroutine la_qconstrained_lstsq_space(A,C,lwork,err)
         !> Least-squares matrix a[m,n]
         real(qp),intent(in) :: A(:,:)
         !> Equality constraint matrix c[p,n]
         real(qp),intent(in) :: C(:,:)
         !> Size of the working space array
         integer(ilp),intent(out) :: lwork
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,p,info
         real(qp) :: a_dummy(1,1),b_dummy(1)
         real(qp) :: c_dummy(1,1),d_dummy(1)
         real(qp) :: work(1),x(1)

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         p = size(C,1,kind=ilp)

         lwork = -1_ilp
         call gglse(m,n,p,a_dummy,m,c_dummy,p,b_dummy,d_dummy,x,work,lwork,info)
         call la_handle_gglse(info,m,n,p,err0)

         lwork = ceiling(real(work(1),kind=qp),kind=ilp)

         call err0%handle(err)

     end subroutine la_qconstrained_lstsq_space

     !> Solve the real(qp) equality-constrained least-squares problem into x
     subroutine la_qsolve_constrained_lstsq(A,b,C,d,x,storage,overwrite_matrices,err)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         real(qp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         real(qp),intent(inout),target :: C(:,:),d(:)
         !> Solution vector x[n]
         real(qp),intent(out) :: x(:)
         !> [optional] working storage space
         real(qp),optional,intent(out),target :: storage(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: ma,na,mb,mc,nc,md,mx,lwork,info
         logical(lk) :: copy_matrices
         real(qp),pointer :: amat(:,:),bvec(:),cmat(:,:),dvec(:),work(:)
         character(*),parameter :: this = 'constrained_lstsq'

         ma = size(A,1,kind=ilp)
         na = size(A,2,kind=ilp)
         mc = size(C,1,kind=ilp)
         nc = size(C,2,kind=ilp)
         mb = size(b,kind=ilp)
         md = size(d,kind=ilp)
         mx = size(x,kind=ilp)

         call la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err0)
         if (err0%error()) then
            call err0%handle(err)
            return
         end if

         ! Can the input matrices be overwritten? By default, do not overwrite
         if (present(overwrite_matrices)) then
            copy_matrices = .not. overwrite_matrices
         else
            copy_matrices = .true._lk
         end if

         if (copy_matrices) then
            allocate (amat(ma,na),source=A)
            allocate (bvec(mb),source=b)
            allocate (cmat(mc,nc),source=C)
            allocate (dvec(md),source=d)
         else
            amat => A
            bvec => b
            cmat => C
            dvec => d
         end if

         call la_qconstrained_lstsq_space(A,C,lwork,err0)

         if (err0%ok()) then

            if (present(storage)) then
               work => storage
            else
               allocate (work(lwork))
            end if

            if (size(work,kind=ilp) < lwork) then
               err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:',size(work,kind=ilp), &
                                                       ' should be >=',lwork)
            else
               call gglse(ma,na,mc,amat,ma,cmat,mc,bvec,dvec,x,work,lwork,info)
               call la_handle_gglse(info,ma,na,mc,err0)
            end if

            if (.not. present(storage)) deallocate (work)

         end if

         if (copy_matrices) deallocate (amat,bvec,cmat,dvec)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_constrained_lstsq

     !> Solve the real(qp) equality-constrained least-squares problem
     function la_qconstrained_lstsq(A,b,C,d,overwrite_matrices,err) result(x)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         real(qp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         real(qp),intent(inout),target :: C(:,:),d(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Solution vector x[n]
         real(qp),allocatable,target :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(A,2,kind=ilp)
         allocate (x(n))

         call la_qsolve_constrained_lstsq(A,b,C,d,x,overwrite_matrices=overwrite_matrices,err=err)

     end function la_qconstrained_lstsq

     !> Compute the real(qp) weighted least-squares solution
     function la_qweighted_lstsq(w,a,b,cond,overwrite_a,rank,err) result(x)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(qp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         real(qp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array x[n]
         real(qp),allocatable :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(a,2,kind=ilp)
         allocate (x(n))

         call la_qsolve_weighted_lstsq(w,a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err)

     end function la_qweighted_lstsq

     !> Compute the real(qp) weighted least-squares solution into x: minimize ||D(Ax - b)||, D = diag(sqrt(w))
     subroutine la_qsolve_weighted_lstsq(w,a,b,x,cond,overwrite_a,rank,err)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(qp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         real(qp),intent(in) :: b(:)
         !> Result array x[n]
         real(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,i
         logical(lk) :: copy_a
         real(qp),pointer :: amat(:,:)
         real(qp),allocatable,target :: amat_alloc(:,:)
         real(qp),allocatable :: b_scaled(:)
         real(qp),allocatable :: sqrt_w(:)
         character(*),parameter :: this = 'weighted_lstsq'

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m < 1 .or. n < 1) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [m,n])
         elseif (size(w,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(w) =',size(w,kind=ilp),'but a has',m,'rows')
         elseif (size(b,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(b) =',size(b,kind=ilp),'but a has',m,'rows')
         elseif (any(w <= 0.0_qp)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'the weights must be positive')
         end if

         if (err0%error()) then
            call err0%handle(err)
            if (present(rank)) rank = 0_ilp
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (copy_a) then
            allocate (amat_alloc(m,n),source=a)
            amat => amat_alloc
         else
            amat => a
         end if

         ! Scale the rows of A and of b by the square root of the weights
         sqrt_w = sqrt(w)
         do i = 1,m
            amat(i,:) = sqrt_w(i)*amat(i,:)
         end do
         b_scaled = sqrt_w*b

         call la_qsolve_lstsq_one(amat,b_scaled,x,cond=cond,overwrite_a=.true._lk,rank=rank,err=err0)

         ! Report the error against this procedure, not against the transformed problem
         if (err0%error()) err0%where_at = this

         if (copy_a) deallocate (amat_alloc)
         deallocate (b_scaled,sqrt_w)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_qsolve_weighted_lstsq

     !> Working space needed by the complex(sp) equality-constrained least-squares solver
     subroutine la_cconstrained_lstsq_space(A,C,lwork,err)
         !> Least-squares matrix a[m,n]
         complex(sp),intent(in) :: A(:,:)
         !> Equality constraint matrix c[p,n]
         complex(sp),intent(in) :: C(:,:)
         !> Size of the working space array
         integer(ilp),intent(out) :: lwork
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,p,info
         complex(sp) :: a_dummy(1,1),b_dummy(1)
         complex(sp) :: c_dummy(1,1),d_dummy(1)
         complex(sp) :: work(1),x(1)

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         p = size(C,1,kind=ilp)

         lwork = -1_ilp
         call gglse(m,n,p,a_dummy,m,c_dummy,p,b_dummy,d_dummy,x,work,lwork,info)
         call la_handle_gglse(info,m,n,p,err0)

         lwork = ceiling(real(work(1),kind=sp),kind=ilp)

         call err0%handle(err)

     end subroutine la_cconstrained_lstsq_space

     !> Solve the complex(sp) equality-constrained least-squares problem into x
     subroutine la_csolve_constrained_lstsq(A,b,C,d,x,storage,overwrite_matrices,err)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         complex(sp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         complex(sp),intent(inout),target :: C(:,:),d(:)
         !> Solution vector x[n]
         complex(sp),intent(out) :: x(:)
         !> [optional] working storage space
         complex(sp),optional,intent(out),target :: storage(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: ma,na,mb,mc,nc,md,mx,lwork,info
         logical(lk) :: copy_matrices
         complex(sp),pointer :: amat(:,:),bvec(:),cmat(:,:),dvec(:),work(:)
         character(*),parameter :: this = 'constrained_lstsq'

         ma = size(A,1,kind=ilp)
         na = size(A,2,kind=ilp)
         mc = size(C,1,kind=ilp)
         nc = size(C,2,kind=ilp)
         mb = size(b,kind=ilp)
         md = size(d,kind=ilp)
         mx = size(x,kind=ilp)

         call la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err0)
         if (err0%error()) then
            call err0%handle(err)
            return
         end if

         ! Can the input matrices be overwritten? By default, do not overwrite
         if (present(overwrite_matrices)) then
            copy_matrices = .not. overwrite_matrices
         else
            copy_matrices = .true._lk
         end if

         if (copy_matrices) then
            allocate (amat(ma,na),source=A)
            allocate (bvec(mb),source=b)
            allocate (cmat(mc,nc),source=C)
            allocate (dvec(md),source=d)
         else
            amat => A
            bvec => b
            cmat => C
            dvec => d
         end if

         call la_cconstrained_lstsq_space(A,C,lwork,err0)

         if (err0%ok()) then

            if (present(storage)) then
               work => storage
            else
               allocate (work(lwork))
            end if

            if (size(work,kind=ilp) < lwork) then
               err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:',size(work,kind=ilp), &
                                                       ' should be >=',lwork)
            else
               call gglse(ma,na,mc,amat,ma,cmat,mc,bvec,dvec,x,work,lwork,info)
               call la_handle_gglse(info,ma,na,mc,err0)
            end if

            if (.not. present(storage)) deallocate (work)

         end if

         if (copy_matrices) deallocate (amat,bvec,cmat,dvec)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_constrained_lstsq

     !> Solve the complex(sp) equality-constrained least-squares problem
     function la_cconstrained_lstsq(A,b,C,d,overwrite_matrices,err) result(x)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         complex(sp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         complex(sp),intent(inout),target :: C(:,:),d(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Solution vector x[n]
         complex(sp),allocatable,target :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(A,2,kind=ilp)
         allocate (x(n))

         call la_csolve_constrained_lstsq(A,b,C,d,x,overwrite_matrices=overwrite_matrices,err=err)

     end function la_cconstrained_lstsq

     !> Compute the complex(sp) weighted least-squares solution
     function la_cweighted_lstsq(w,a,b,cond,overwrite_a,rank,err) result(x)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(sp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         complex(sp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array x[n]
         complex(sp),allocatable :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(a,2,kind=ilp)
         allocate (x(n))

         call la_csolve_weighted_lstsq(w,a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err)

     end function la_cweighted_lstsq

     !> Compute the complex(sp) weighted least-squares solution into x: minimize ||D(Ax - b)||, D = diag(sqrt(w))
     subroutine la_csolve_weighted_lstsq(w,a,b,x,cond,overwrite_a,rank,err)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(sp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         complex(sp),intent(in) :: b(:)
         !> Result array x[n]
         complex(sp),intent(inout),contiguous,target :: x(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(sp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,i
         logical(lk) :: copy_a
         complex(sp),pointer :: amat(:,:)
         complex(sp),allocatable,target :: amat_alloc(:,:)
         complex(sp),allocatable :: b_scaled(:)
         real(sp),allocatable :: sqrt_w(:)
         character(*),parameter :: this = 'weighted_lstsq'

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m < 1 .or. n < 1) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [m,n])
         elseif (size(w,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(w) =',size(w,kind=ilp),'but a has',m,'rows')
         elseif (size(b,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(b) =',size(b,kind=ilp),'but a has',m,'rows')
         elseif (any(w <= 0.0_sp)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'the weights must be positive')
         end if

         if (err0%error()) then
            call err0%handle(err)
            if (present(rank)) rank = 0_ilp
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (copy_a) then
            allocate (amat_alloc(m,n),source=a)
            amat => amat_alloc
         else
            amat => a
         end if

         ! Scale the rows of A and of b by the square root of the weights
         sqrt_w = sqrt(w)
         do i = 1,m
            amat(i,:) = sqrt_w(i)*amat(i,:)
         end do
         b_scaled = sqrt_w*b

         call la_csolve_lstsq_one(amat,b_scaled,x,cond=cond,overwrite_a=.true._lk,rank=rank,err=err0)

         ! Report the error against this procedure, not against the transformed problem
         if (err0%error()) err0%where_at = this

         if (copy_a) deallocate (amat_alloc)
         deallocate (b_scaled,sqrt_w)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_csolve_weighted_lstsq

     !> Working space needed by the complex(dp) equality-constrained least-squares solver
     subroutine la_zconstrained_lstsq_space(A,C,lwork,err)
         !> Least-squares matrix a[m,n]
         complex(dp),intent(in) :: A(:,:)
         !> Equality constraint matrix c[p,n]
         complex(dp),intent(in) :: C(:,:)
         !> Size of the working space array
         integer(ilp),intent(out) :: lwork
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,p,info
         complex(dp) :: a_dummy(1,1),b_dummy(1)
         complex(dp) :: c_dummy(1,1),d_dummy(1)
         complex(dp) :: work(1),x(1)

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         p = size(C,1,kind=ilp)

         lwork = -1_ilp
         call gglse(m,n,p,a_dummy,m,c_dummy,p,b_dummy,d_dummy,x,work,lwork,info)
         call la_handle_gglse(info,m,n,p,err0)

         lwork = ceiling(real(work(1),kind=dp),kind=ilp)

         call err0%handle(err)

     end subroutine la_zconstrained_lstsq_space

     !> Solve the complex(dp) equality-constrained least-squares problem into x
     subroutine la_zsolve_constrained_lstsq(A,b,C,d,x,storage,overwrite_matrices,err)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         complex(dp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         complex(dp),intent(inout),target :: C(:,:),d(:)
         !> Solution vector x[n]
         complex(dp),intent(out) :: x(:)
         !> [optional] working storage space
         complex(dp),optional,intent(out),target :: storage(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: ma,na,mb,mc,nc,md,mx,lwork,info
         logical(lk) :: copy_matrices
         complex(dp),pointer :: amat(:,:),bvec(:),cmat(:,:),dvec(:),work(:)
         character(*),parameter :: this = 'constrained_lstsq'

         ma = size(A,1,kind=ilp)
         na = size(A,2,kind=ilp)
         mc = size(C,1,kind=ilp)
         nc = size(C,2,kind=ilp)
         mb = size(b,kind=ilp)
         md = size(d,kind=ilp)
         mx = size(x,kind=ilp)

         call la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err0)
         if (err0%error()) then
            call err0%handle(err)
            return
         end if

         ! Can the input matrices be overwritten? By default, do not overwrite
         if (present(overwrite_matrices)) then
            copy_matrices = .not. overwrite_matrices
         else
            copy_matrices = .true._lk
         end if

         if (copy_matrices) then
            allocate (amat(ma,na),source=A)
            allocate (bvec(mb),source=b)
            allocate (cmat(mc,nc),source=C)
            allocate (dvec(md),source=d)
         else
            amat => A
            bvec => b
            cmat => C
            dvec => d
         end if

         call la_zconstrained_lstsq_space(A,C,lwork,err0)

         if (err0%ok()) then

            if (present(storage)) then
               work => storage
            else
               allocate (work(lwork))
            end if

            if (size(work,kind=ilp) < lwork) then
               err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:',size(work,kind=ilp), &
                                                       ' should be >=',lwork)
            else
               call gglse(ma,na,mc,amat,ma,cmat,mc,bvec,dvec,x,work,lwork,info)
               call la_handle_gglse(info,ma,na,mc,err0)
            end if

            if (.not. present(storage)) deallocate (work)

         end if

         if (copy_matrices) deallocate (amat,bvec,cmat,dvec)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_constrained_lstsq

     !> Solve the complex(dp) equality-constrained least-squares problem
     function la_zconstrained_lstsq(A,b,C,d,overwrite_matrices,err) result(x)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         complex(dp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         complex(dp),intent(inout),target :: C(:,:),d(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Solution vector x[n]
         complex(dp),allocatable,target :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(A,2,kind=ilp)
         allocate (x(n))

         call la_zsolve_constrained_lstsq(A,b,C,d,x,overwrite_matrices=overwrite_matrices,err=err)

     end function la_zconstrained_lstsq

     !> Compute the complex(dp) weighted least-squares solution
     function la_zweighted_lstsq(w,a,b,cond,overwrite_a,rank,err) result(x)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(dp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         complex(dp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array x[n]
         complex(dp),allocatable :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(a,2,kind=ilp)
         allocate (x(n))

         call la_zsolve_weighted_lstsq(w,a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err)

     end function la_zweighted_lstsq

     !> Compute the complex(dp) weighted least-squares solution into x: minimize ||D(Ax - b)||, D = diag(sqrt(w))
     subroutine la_zsolve_weighted_lstsq(w,a,b,x,cond,overwrite_a,rank,err)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(dp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         complex(dp),intent(in) :: b(:)
         !> Result array x[n]
         complex(dp),intent(inout),contiguous,target :: x(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(dp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,i
         logical(lk) :: copy_a
         complex(dp),pointer :: amat(:,:)
         complex(dp),allocatable,target :: amat_alloc(:,:)
         complex(dp),allocatable :: b_scaled(:)
         real(dp),allocatable :: sqrt_w(:)
         character(*),parameter :: this = 'weighted_lstsq'

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m < 1 .or. n < 1) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [m,n])
         elseif (size(w,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(w) =',size(w,kind=ilp),'but a has',m,'rows')
         elseif (size(b,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(b) =',size(b,kind=ilp),'but a has',m,'rows')
         elseif (any(w <= 0.0_dp)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'the weights must be positive')
         end if

         if (err0%error()) then
            call err0%handle(err)
            if (present(rank)) rank = 0_ilp
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (copy_a) then
            allocate (amat_alloc(m,n),source=a)
            amat => amat_alloc
         else
            amat => a
         end if

         ! Scale the rows of A and of b by the square root of the weights
         sqrt_w = sqrt(w)
         do i = 1,m
            amat(i,:) = sqrt_w(i)*amat(i,:)
         end do
         b_scaled = sqrt_w*b

         call la_zsolve_lstsq_one(amat,b_scaled,x,cond=cond,overwrite_a=.true._lk,rank=rank,err=err0)

         ! Report the error against this procedure, not against the transformed problem
         if (err0%error()) err0%where_at = this

         if (copy_a) deallocate (amat_alloc)
         deallocate (b_scaled,sqrt_w)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_zsolve_weighted_lstsq

     !> Working space needed by the complex(qp) equality-constrained least-squares solver
     subroutine la_wconstrained_lstsq_space(A,C,lwork,err)
         !> Least-squares matrix a[m,n]
         complex(qp),intent(in) :: A(:,:)
         !> Equality constraint matrix c[p,n]
         complex(qp),intent(in) :: C(:,:)
         !> Size of the working space array
         integer(ilp),intent(out) :: lwork
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,p,info
         complex(qp) :: a_dummy(1,1),b_dummy(1)
         complex(qp) :: c_dummy(1,1),d_dummy(1)
         complex(qp) :: work(1),x(1)

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         p = size(C,1,kind=ilp)

         lwork = -1_ilp
         call gglse(m,n,p,a_dummy,m,c_dummy,p,b_dummy,d_dummy,x,work,lwork,info)
         call la_handle_gglse(info,m,n,p,err0)

         lwork = ceiling(real(work(1),kind=qp),kind=ilp)

         call err0%handle(err)

     end subroutine la_wconstrained_lstsq_space

     !> Solve the complex(qp) equality-constrained least-squares problem into x
     subroutine la_wsolve_constrained_lstsq(A,b,C,d,x,storage,overwrite_matrices,err)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         complex(qp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         complex(qp),intent(inout),target :: C(:,:),d(:)
         !> Solution vector x[n]
         complex(qp),intent(out) :: x(:)
         !> [optional] working storage space
         complex(qp),optional,intent(out),target :: storage(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: ma,na,mb,mc,nc,md,mx,lwork,info
         logical(lk) :: copy_matrices
         complex(qp),pointer :: amat(:,:),bvec(:),cmat(:,:),dvec(:),work(:)
         character(*),parameter :: this = 'constrained_lstsq'

         ma = size(A,1,kind=ilp)
         na = size(A,2,kind=ilp)
         mc = size(C,1,kind=ilp)
         nc = size(C,2,kind=ilp)
         mb = size(b,kind=ilp)
         md = size(d,kind=ilp)
         mx = size(x,kind=ilp)

         call la_check_constrained_sizes(ma,na,mb,mc,nc,md,mx,err0)
         if (err0%error()) then
            call err0%handle(err)
            return
         end if

         ! Can the input matrices be overwritten? By default, do not overwrite
         if (present(overwrite_matrices)) then
            copy_matrices = .not. overwrite_matrices
         else
            copy_matrices = .true._lk
         end if

         if (copy_matrices) then
            allocate (amat(ma,na),source=A)
            allocate (bvec(mb),source=b)
            allocate (cmat(mc,nc),source=C)
            allocate (dvec(md),source=d)
         else
            amat => A
            bvec => b
            cmat => C
            dvec => d
         end if

         call la_wconstrained_lstsq_space(A,C,lwork,err0)

         if (err0%ok()) then

            if (present(storage)) then
               work => storage
            else
               allocate (work(lwork))
            end if

            if (size(work,kind=ilp) < lwork) then
               err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient working space:',size(work,kind=ilp), &
                                                       ' should be >=',lwork)
            else
               call gglse(ma,na,mc,amat,ma,cmat,mc,bvec,dvec,x,work,lwork,info)
               call la_handle_gglse(info,ma,na,mc,err0)
            end if

            if (.not. present(storage)) deallocate (work)

         end if

         if (copy_matrices) deallocate (amat,bvec,cmat,dvec)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_constrained_lstsq

     !> Solve the complex(qp) equality-constrained least-squares problem
     function la_wconstrained_lstsq(A,b,C,d,overwrite_matrices,err) result(x)
         !> Least-squares matrix a[m,n] and right hand side b[m]
         complex(qp),intent(inout),target :: A(:,:),b(:)
         !> Equality constraint matrix c[p,n] and right hand side d[p]
         complex(qp),intent(inout),target :: C(:,:),d(:)
         !> [optional] Can A,b,C,d data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_matrices
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Solution vector x[n]
         complex(qp),allocatable,target :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(A,2,kind=ilp)
         allocate (x(n))

         call la_wsolve_constrained_lstsq(A,b,C,d,x,overwrite_matrices=overwrite_matrices,err=err)

     end function la_wconstrained_lstsq

     !> Compute the complex(qp) weighted least-squares solution
     function la_wweighted_lstsq(w,a,b,cond,overwrite_a,rank,err) result(x)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(qp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         complex(qp),intent(in) :: b(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         !> Result array x[n]
         complex(qp),allocatable :: x(:)

         !> Local variables
         integer(ilp) :: n

         n = size(a,2,kind=ilp)
         allocate (x(n))

         call la_wsolve_weighted_lstsq(w,a,b,x,cond=cond,overwrite_a=overwrite_a,rank=rank,err=err)

     end function la_wweighted_lstsq

     !> Compute the complex(qp) weighted least-squares solution into x: minimize ||D(Ax - b)||, D = diag(sqrt(w))
     subroutine la_wsolve_weighted_lstsq(w,a,b,x,cond,overwrite_a,rank,err)
         !> Weight vector w[m]. It is always real, and all its entries must be positive
         real(qp),intent(in) :: w(:)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Right hand side vector b[m]
         complex(qp),intent(in) :: b(:)
         !> Result array x[n]
         complex(qp),intent(inout),contiguous,target :: x(:)
         !> [optional] cutoff for rank evaluation: singular values s(i)<=cond*maxval(s) are considered 0.
         real(qp),optional,intent(in) :: cond
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] Return rank of A
         integer(ilp),optional,intent(out) :: rank
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: m,n,i
         logical(lk) :: copy_a
         complex(qp),pointer :: amat(:,:)
         complex(qp),allocatable,target :: amat_alloc(:,:)
         complex(qp),allocatable :: b_scaled(:)
         real(qp),allocatable :: sqrt_w(:)
         character(*),parameter :: this = 'weighted_lstsq'

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         if (m < 1 .or. n < 1) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a(m,n) =', [m,n])
         elseif (size(w,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(w) =',size(w,kind=ilp),'but a has',m,'rows')
         elseif (size(b,kind=ilp) /= m) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'size(b) =',size(b,kind=ilp),'but a has',m,'rows')
         elseif (any(w <= 0.0_qp)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'the weights must be positive')
         end if

         if (err0%error()) then
            call err0%handle(err)
            if (present(rank)) rank = 0_ilp
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            copy_a = .not. overwrite_a
         else
            copy_a = .true._lk
         end if

         if (copy_a) then
            allocate (amat_alloc(m,n),source=a)
            amat => amat_alloc
         else
            amat => a
         end if

         ! Scale the rows of A and of b by the square root of the weights
         sqrt_w = sqrt(w)
         do i = 1,m
            amat(i,:) = sqrt_w(i)*amat(i,:)
         end do
         b_scaled = sqrt_w*b

         call la_wsolve_lstsq_one(amat,b_scaled,x,cond=cond,overwrite_a=.true._lk,rank=rank,err=err0)

         ! Report the error against this procedure, not against the transformed problem
         if (err0%error()) err0%where_at = this

         if (copy_a) deallocate (amat_alloc)
         deallocate (b_scaled,sqrt_w)

         ! Process output and return
         call err0%handle(err)

     end subroutine la_wsolve_weighted_lstsq

end module la_least_squares
