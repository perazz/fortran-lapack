!> Functions of a square matrix
module la_matrix_functions
     use la_constants
     use la_blas,only:gemm
     use la_lapack,only:gesv
     use la_state_type
     implicit none(type,external)
     private

     public :: expm
     public :: matrix_exp
     public :: expm_space

     character(*),parameter :: this = 'matrix_exponential'

     !> @brief Compute the exponential of a square matrix.
     !!
     !! This function computes \f$ E = e^{A} \f$ for a real or complex square matrix \f$ A \f$,
     !! by scaling and squaring with a Pade approximant.
     !!
     !! @param[in] A The input square matrix of size \f$ [n, n] \f$.
     !! @param[in] order (Optional) Order of the Pade approximant. Default is 10.
     !! @param[out] err (Optional) A state return flag. Without it the function is `pure`, and an
     !!                 error stops execution.
     !!
     !! @return The matrix exponential \f$ e^{A} \f$, of the same kind and size as \f$ A \f$.
     !!
     !! @note The function form without `err` is `pure`; the form with `err` is not, because a
     !!       `pure` function cannot return an error argument.
     !!
     interface expm
        module procedure la_s_expm_fun
        module procedure la_s_expm_fun_err
        module procedure la_d_expm_fun
        module procedure la_d_expm_fun_err
#ifdef LA_WITH_XDP
        module procedure la_x_expm_fun
        module procedure la_x_expm_fun_err
#endif
#ifdef LA_WITH_QP
        module procedure la_q_expm_fun
        module procedure la_q_expm_fun_err
#endif
        module procedure la_c_expm_fun
        module procedure la_c_expm_fun_err
        module procedure la_z_expm_fun
        module procedure la_z_expm_fun_err
#ifdef LA_WITH_XDP
        module procedure la_y_expm_fun
        module procedure la_y_expm_fun_err
#endif
#ifdef LA_WITH_QP
        module procedure la_w_expm_fun
        module procedure la_w_expm_fun_err
#endif
     end interface expm

     !> @brief Compute the exponential of a square matrix in a subroutine.
     !!
     !! This subroutine computes \f$ E = e^{A} \f$ for a real or complex square matrix \f$ A \f$,
     !! by scaling and squaring with a Pade approximant. It either overwrites \f$ A \f$ with its
     !! exponential, or writes the exponential into a second matrix \f$ E \f$.
     !!
     !! @param[in,out] A The input square matrix of size \f$ [n, n] \f$. Without `E`, it is
     !!                  overwritten with \f$ e^{A} \f$.
     !! @param[out] E (Optional) The output matrix of size \f$ [n, n] \f$.
     !! @param[in] order (Optional) Order of the Pade approximant. Default is 10.
     !! @param[in,out] work (Optional) Workspace of the same type and kind as \f$ A \f$, of size at
     !!                     least the `lwork` returned by [expm_space](@ref la_matrix_functions::expm_space).
     !! @param[in,out] iwork (Optional) Integer workspace of size at least the `liwork` returned by
     !!                      [expm_space](@ref la_matrix_functions::expm_space).
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not
     !!                 provided, the subroutine will stop execution.
     !!
     !! @note The subroutine is `pure`. With both `work` and `iwork` present it allocates no work
     !!       arrays of its own.
     !!
     interface matrix_exp
        module procedure la_s_expm_inplace
        module procedure la_s_expm
        module procedure la_d_expm_inplace
        module procedure la_d_expm
#ifdef LA_WITH_XDP
        module procedure la_x_expm_inplace
        module procedure la_x_expm
#endif
#ifdef LA_WITH_QP
        module procedure la_q_expm_inplace
        module procedure la_q_expm
#endif
        module procedure la_c_expm_inplace
        module procedure la_c_expm
        module procedure la_z_expm_inplace
        module procedure la_z_expm
#ifdef LA_WITH_XDP
        module procedure la_y_expm_inplace
        module procedure la_y_expm
#endif
#ifdef LA_WITH_QP
        module procedure la_w_expm_inplace
        module procedure la_w_expm
#endif
     end interface matrix_exp

     !> @brief Workspace sizes for the matrix exponential.
     !!
     !! This subroutine returns the sizes of the working arrays that
     !! [matrix_exp](@ref la_matrix_functions::matrix_exp) needs for a matrix of the size of \f$ A \f$,
     !! so that repeated calls perform no internal allocation.
     !!
     !! @param[in] A The input square matrix of size \f$ [n, n] \f$. Only its size is read.
     !! @param[out] lwork The minimum size of `work`, \f$ 4 n^2 \f$.
     !! @param[out] liwork The minimum size of `iwork`, \f$ n \f$.
     !!
     interface expm_space
        module procedure la_s_expm_space
        module procedure la_d_expm_space
#ifdef LA_WITH_XDP
        module procedure la_x_expm_space
#endif
#ifdef LA_WITH_QP
        module procedure la_q_expm_space
#endif
        module procedure la_c_expm_space
        module procedure la_z_expm_space
#ifdef LA_WITH_XDP
        module procedure la_y_expm_space
#endif
#ifdef LA_WITH_QP
        module procedure la_w_expm_space
#endif
     end interface expm_space

     contains

     elemental subroutine handle_gesv_info(info,n,err)
         integer(ilp),intent(in) :: info,n
         type(la_state),intent(out) :: err

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (1:)
                err = la_state(this,LINALG_ERROR,'singular Pade denominator')
            case default
                err = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gesv_info

     !> Work array sizes for a matrix of size n: four n*n matrices and n pivots
     pure subroutine expm_workspace(n,lwork,liwork)
         integer(ilp),intent(in) :: n
         integer(ilp),intent(out) :: lwork,liwork

         lwork = 4*n*n
         liwork = n

     end subroutine expm_workspace

     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_s_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         real(sp),allocatable :: e(:,:)

         e = a
         call la_s_expm_inplace(e,order)

     end function la_s_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_s_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         real(sp),allocatable :: e(:,:)

         e = a
         call la_s_expm_inplace(e,order,err=err)

     end function la_s_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_s_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         real(sp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_s_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_s_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         real(sp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(sp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_s_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_s_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_s_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         real(sp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(sp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         real(sp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_s_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_s_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_s_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         real(sp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         real(sp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         real(sp),parameter :: zero = 0.0_sp
         real(sp),parameter :: one = 1.0_sp
         integer(ilp) :: ee,k,s,i,j
         real(sp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_sp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_sp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_sp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_sp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_sp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_sp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_sp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_s_expm_kernel

     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_d_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         real(dp),allocatable :: e(:,:)

         e = a
         call la_d_expm_inplace(e,order)

     end function la_d_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_d_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         real(dp),allocatable :: e(:,:)

         e = a
         call la_d_expm_inplace(e,order,err=err)

     end function la_d_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_d_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         real(dp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_d_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_d_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         real(dp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(dp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_d_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_d_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_d_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         real(dp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(dp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         real(dp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_d_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_d_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_d_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         real(dp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         real(dp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         real(dp),parameter :: zero = 0.0_dp
         real(dp),parameter :: one = 1.0_dp
         integer(ilp) :: ee,k,s,i,j
         real(dp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_dp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_dp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_dp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_dp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_dp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_dp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_dp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_d_expm_kernel

#ifdef LA_WITH_XDP
     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_x_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         real(xdp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         real(xdp),allocatable :: e(:,:)

         e = a
         call la_x_expm_inplace(e,order)

     end function la_x_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_x_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         real(xdp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         real(xdp),allocatable :: e(:,:)

         e = a
         call la_x_expm_inplace(e,order,err=err)

     end function la_x_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_x_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         real(xdp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_x_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_x_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         real(xdp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         real(xdp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(xdp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_x_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_x_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_x_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         real(xdp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(xdp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         real(xdp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_x_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_x_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_x_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         real(xdp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         real(xdp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         real(xdp),parameter :: zero = 0.0_xdp
         real(xdp),parameter :: one = 1.0_xdp
         integer(ilp) :: ee,k,s,i,j
         real(xdp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_xdp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_xdp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_xdp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_xdp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_xdp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_xdp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_xdp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_x_expm_kernel
#endif

#ifdef LA_WITH_QP
     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_q_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         real(qp),allocatable :: e(:,:)

         e = a
         call la_q_expm_inplace(e,order)

     end function la_q_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_q_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         real(qp),allocatable :: e(:,:)

         e = a
         call la_q_expm_inplace(e,order,err=err)

     end function la_q_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_q_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         real(qp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_q_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_q_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         real(qp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(qp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_q_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_q_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_q_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         real(qp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         real(qp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         real(qp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_q_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_q_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_q_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         real(qp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         real(qp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         real(qp),parameter :: zero = 0.0_qp
         real(qp),parameter :: one = 1.0_qp
         integer(ilp) :: ee,k,s,i,j
         real(qp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_qp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_qp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_qp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_qp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_qp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_qp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_qp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_q_expm_kernel
#endif

     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_c_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         complex(sp),allocatable :: e(:,:)

         e = a
         call la_c_expm_inplace(e,order)

     end function la_c_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_c_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         complex(sp),allocatable :: e(:,:)

         e = a
         call la_c_expm_inplace(e,order,err=err)

     end function la_c_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_c_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         complex(sp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_c_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_c_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         complex(sp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(sp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_c_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_c_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_c_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         complex(sp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(sp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         complex(sp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_c_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_c_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_c_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         complex(sp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         complex(sp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         complex(sp),parameter :: zero = 0.0_sp
         complex(sp),parameter :: one = 1.0_sp
         integer(ilp) :: ee,k,s,i,j
         real(sp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_sp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_sp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_sp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_sp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_sp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_sp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_sp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_c_expm_kernel

     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_z_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         complex(dp),allocatable :: e(:,:)

         e = a
         call la_z_expm_inplace(e,order)

     end function la_z_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_z_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         complex(dp),allocatable :: e(:,:)

         e = a
         call la_z_expm_inplace(e,order,err=err)

     end function la_z_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_z_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         complex(dp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_z_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_z_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         complex(dp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(dp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_z_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_z_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_z_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         complex(dp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(dp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         complex(dp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_z_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_z_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_z_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         complex(dp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         complex(dp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         complex(dp),parameter :: zero = 0.0_dp
         complex(dp),parameter :: one = 1.0_dp
         integer(ilp) :: ee,k,s,i,j
         real(dp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_dp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_dp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_dp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_dp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_dp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_dp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_dp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_z_expm_kernel

#ifdef LA_WITH_XDP
     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_y_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         complex(xdp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         complex(xdp),allocatable :: e(:,:)

         e = a
         call la_y_expm_inplace(e,order)

     end function la_y_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_y_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         complex(xdp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         complex(xdp),allocatable :: e(:,:)

         e = a
         call la_y_expm_inplace(e,order,err=err)

     end function la_y_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_y_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         complex(xdp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_y_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_y_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         complex(xdp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         complex(xdp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(xdp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_y_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_y_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_y_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         complex(xdp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(xdp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         complex(xdp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_y_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_y_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_y_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         complex(xdp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         complex(xdp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         complex(xdp),parameter :: zero = 0.0_xdp
         complex(xdp),parameter :: one = 1.0_xdp
         integer(ilp) :: ee,k,s,i,j
         real(xdp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_xdp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_xdp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_xdp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_xdp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_xdp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_xdp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_xdp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_y_expm_kernel
#endif

#ifdef LA_WITH_QP
     !> Matrix exponential, function form: `pure`, stops on error
     pure function la_w_expm_fun(a,order) result(e)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> Matrix exponential e[n,n]
         complex(qp),allocatable :: e(:,:)

         e = a
         call la_w_expm_inplace(e,order)

     end function la_w_expm_fun

     !> Matrix exponential, function form with a state return flag
     function la_w_expm_fun_err(a,order,err) result(e)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> State return flag
         type(la_state),intent(out) :: err
         !> Matrix exponential e[n,n]
         complex(qp),allocatable :: e(:,:)

         e = a
         call la_w_expm_inplace(e,order,err=err)

     end function la_w_expm_fun_err

     !> Workspace sizes for the matrix exponential of a[n,n]
     pure subroutine la_w_expm_space(a,lwork,liwork)
         !> Input matrix a[n,n]; only its size is read
         complex(qp),intent(in) :: a(:,:)
         !> Minimum size of the work array
         integer(ilp),intent(out) :: lwork
         !> Minimum size of the integer work array
         integer(ilp),intent(out) :: liwork

         call expm_workspace(size(a,2,kind=ilp),lwork,liwork)

     end subroutine la_w_expm_space

     !> Matrix exponential written into a second matrix
     pure subroutine la_w_expm(a,e,order,work,iwork,err)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Output matrix exponential e[n,n]
         complex(qp),intent(out) :: e(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(qp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: lda,n,lde,ne

         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         lde = size(e,1,kind=ilp)
         ne = size(e,2,kind=ilp)

         if (lda /= n .or. lde /= n .or. ne /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],', e=', [lde,ne])
         else
            e = a
            call la_w_expm_inplace(e,order,work,iwork,err0)
         end if

         call err0%handle(err)

     end subroutine la_w_expm

     !> Matrix exponential in place: a is overwritten with exp(a)
     pure subroutine la_w_expm_inplace(a,order,work,iwork,err)
         !> Input matrix a[n,n]; on output, its exponential
         complex(qp),intent(inout) :: a(:,:)
         !> [optional] Order of the Pade approximant
         integer(ilp),optional,intent(in) :: order
         !> [optional] Workspace, size to be checked with expm_space
         complex(qp),optional,intent(inout),contiguous,target :: work(:)
         !> [optional] Integer workspace, size to be checked with expm_space
         integer(ilp),optional,intent(inout),contiguous,target :: iwork(:)
         !> [optional] State return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         type(la_state) :: err0
         integer(ilp) :: m,n,order_,lwork,liwork,nwork,niwork,info
         integer(ilp),pointer,contiguous :: ipiv(:)
         complex(qp),pointer,contiguous :: wk(:,:,:)

         order_ = 10
         if (present(order)) order_ = order

         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Sizes needed, and sizes of the caller's arrays if any
         call expm_workspace(n,lwork,liwork)
         nwork = lwork
         niwork = liwork
         if (present(work)) nwork = size(work,kind=ilp)
         if (present(iwork)) niwork = size(iwork,kind=ilp)

         if (m /= n) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
         else if (order_ < 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                            'the order of the Pade approximant cannot be negative, order=',order_)
         else if (nwork < lwork .or. niwork < liwork) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'insufficient workspace: work=',nwork, &
                            ' should be >=',lwork,', iwork=',niwork,' should be >=',liwork)
         else if (n > 0) then
            ! User or self-allocated workspace: four n*n matrices and n pivots
            if (present(work)) then
               wk(1:n,1:n,1:4) => work(1:lwork)
            else
               allocate (wk(n,n,4))
            end if
            if (present(iwork)) then
               ipiv => iwork
            else
               allocate (ipiv(liwork))
            end if

            call la_w_expm_kernel(n,order_,a,wk(:,:,1),wk(:,:,2),wk(:,:,3),wk(:,:,4),ipiv,info)
            call handle_gesv_info(info,n,err0)

            if (.not. present(work)) deallocate (wk)
            if (.not. present(iwork)) deallocate (ipiv)
         end if

         call err0%handle(err)

     end subroutine la_w_expm_inplace

     !> Scaling and squaring with a Pade approximant, on explicit-shape work arrays
     pure subroutine la_w_expm_kernel(n,order,a,a2,q,x,x_tmp,ipiv,info)
         !> Matrix size and order of the Pade approximant
         integer(ilp),intent(in) :: n,order
         !> Input matrix a[n,n]; on output, its exponential
         complex(qp),intent(inout) :: a(n,n)
         !> Work arrays: the scaled matrix, the denominator, and the powers of the scaled matrix
         complex(qp),intent(out) :: a2(n,n),q(n,n),x(n,n),x_tmp(n,n)
         !> Pivots of the Pade solve
         integer(ilp),intent(out) :: ipiv(n)
         !> gesv exit code
         integer(ilp),intent(out) :: info

         complex(qp),parameter :: zero = 0.0_qp
         complex(qp),parameter :: one = 1.0_qp
         integer(ilp) :: ee,k,s,i,j
         real(qp) :: a_norm,c
         logical(lk) :: p

         ! Infinity norm, computed here because the LAPACK lange is not pure
         a_norm = 0.0_qp
         do i = 1,n
            a_norm = max(a_norm,sum(abs(a(i,:))))
         end do

         ! Scaling factor: the scaled matrix has a norm below one half. A zero or non-finite
         ! matrix is not scaled.
         if (a_norm > 0.0_qp .and. a_norm <= huge(a_norm)) then
            ee = int(log(a_norm)/log(2.0_qp),kind=ilp) + 1
            s = max(0_ilp,ee + 1)
         else
            s = 0
         end if

         ! Scale the input matrix and initialize the polynomials
         a2 = a/2.0_qp**s
         x = a2

         ! First step of the Pade approximant
         c = 0.5_qp
         do concurrent(i=1:n,j=1:n)
            a(i,j) = merge(1.0_qp + c*a2(i,j),c*a2(i,j),i == j)
            q(i,j) = merge(1.0_qp - c*a2(i,j),-c*a2(i,j),i == j)
         end do

         ! Remaining terms of numerator (a) and denominator (q)
         p = .true.
         do k = 2,order
            c = c*(order - k + 1)/(k*(2*order - k + 1))
            x_tmp = x
            call gemm('N','N',n,n,n,one,a2,n,x_tmp,n,zero,x,n)
            do concurrent(i=1:n,j=1:n)
               a(i,j) = a(i,j) + c*x(i,j)
               q(i,j) = merge(q(i,j) + c*x(i,j),q(i,j) - c*x(i,j),p)
            end do
            p = .not. p
         end do

         ! exp(a2) = inv(q) * a
         call gesv(n,n,q,n,ipiv,a,n,info)

         ! Undo the scaling by repeated squaring
         if (info == 0) then
            do k = 1,s
               x = a
               call gemm('N','N',n,n,n,one,x,n,x,n,zero,a,n)
            end do
         end if

     end subroutine la_w_expm_kernel
#endif

end module la_matrix_functions
