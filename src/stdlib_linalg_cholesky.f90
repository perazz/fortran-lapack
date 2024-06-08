! Cholesky factorization of a matrix, based on LAPACK *POTRF functions
module stdlib_linalg_cholesky
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack,only:potrf
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Cholesky factorization of a matrix
     public :: cholesky

     ! Returns Lower Cholesky factor
     ! NumPy: L = cholesky(a)
     
     ! Returns the Upper or Lower Cholesky factor of A
     ! SciPy: C = cholesky(a, lower=False, overwrite_a=False, check_finite=True)
     
     ! Returns the Upper of Lower Cholesky factor of A, *without setting zeros in the unused part*
     ! SciPy: [C, lower] = cho_factor(a, lower=False, overwrite_a=False, check_finite=True)
     ! SciPy: cho_solve(c_and_lower, b, overwrite_b=False, check_finite=True)
     
     interface cholesky
        module procedure stdlib_linalg_s_cholesky_inplace
        module procedure stdlib_linalg_s_cholesky
        module procedure stdlib_linalg_d_cholesky_inplace
        module procedure stdlib_linalg_d_cholesky
        module procedure stdlib_linalg_q_cholesky_inplace
        module procedure stdlib_linalg_q_cholesky
        module procedure stdlib_linalg_c_cholesky_inplace
        module procedure stdlib_linalg_c_cholesky
        module procedure stdlib_linalg_z_cholesky_inplace
        module procedure stdlib_linalg_z_cholesky
        module procedure stdlib_linalg_w_cholesky_inplace
        module procedure stdlib_linalg_w_cholesky
     end interface cholesky
     
     character(*),parameter :: this = 'cholesky'

     contains

     elemental subroutine handle_potrf_info(info,triangle,lda,n,err)
         character,intent(in) :: triangle
         integer(ilp),intent(in) :: info,lda,n
         type(linalg_state),intent(out) :: err

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid triangle selection: ', &
                                   triangle,'. should be U/L')
            case (-2)
                err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-4)
                err = linalg_state(this,LINALG_VALUE_ERROR,'invalid lda=',lda,': is < n = ',n)
            case (1:)
                err = linalg_state(this,LINALG_ERROR,'cannot complete factorization:',info, &
                                   '-th order leading minor is not positive definite')
            case default
                err = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_potrf_info
     
     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine stdlib_linalg_s_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,i,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         real(sp),parameter :: zero = 0.0_sp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         lower_ = .true.
         if (present(lower)) lower_ = lower
         triangle = merge('L','U',lower_)
         
         !> Check if the unused half of the return matrix should be zeroed out (default).
         !> Otherwise it is unused and will contain garbage.
         other_zeroed_ = .true.
         if (present(other_zeroed)) other_zeroed_ = other_zeroed

         !> Problem size
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Check sizes
         if (n < 1 .or. lda < 1 .or. lda < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a(m,n)=', [lda,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
         ! Compute factorization
         call potrf(triangle,n,a,lda,info)
         call handle_potrf_info(info,triangle,lda,n,err0)
                  
         ! Zero-out the unused part of matrix A
         clean_unused: if (other_zeroed_ .and. err0%ok()) then
             if (lower_) then
                forall (j=2:n) a(:j - 1,j) = zero
             else
                forall (j=1:n - 1) a(j + 1:,j) = zero
             end if
         end if clean_unused

         ! Process output and return
         call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_s_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine stdlib_linalg_s_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(sp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call stdlib_linalg_s_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call linalg_error_handling(err0,err)
         
     end subroutine stdlib_linalg_s_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine stdlib_linalg_d_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,i,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         real(dp),parameter :: zero = 0.0_dp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         lower_ = .true.
         if (present(lower)) lower_ = lower
         triangle = merge('L','U',lower_)
         
         !> Check if the unused half of the return matrix should be zeroed out (default).
         !> Otherwise it is unused and will contain garbage.
         other_zeroed_ = .true.
         if (present(other_zeroed)) other_zeroed_ = other_zeroed

         !> Problem size
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Check sizes
         if (n < 1 .or. lda < 1 .or. lda < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a(m,n)=', [lda,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
         ! Compute factorization
         call potrf(triangle,n,a,lda,info)
         call handle_potrf_info(info,triangle,lda,n,err0)
                  
         ! Zero-out the unused part of matrix A
         clean_unused: if (other_zeroed_ .and. err0%ok()) then
             if (lower_) then
                forall (j=2:n) a(:j - 1,j) = zero
             else
                forall (j=1:n - 1) a(j + 1:,j) = zero
             end if
         end if clean_unused

         ! Process output and return
         call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_d_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine stdlib_linalg_d_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(dp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call stdlib_linalg_d_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call linalg_error_handling(err0,err)
         
     end subroutine stdlib_linalg_d_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine stdlib_linalg_q_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,i,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         real(qp),parameter :: zero = 0.0_qp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         lower_ = .true.
         if (present(lower)) lower_ = lower
         triangle = merge('L','U',lower_)
         
         !> Check if the unused half of the return matrix should be zeroed out (default).
         !> Otherwise it is unused and will contain garbage.
         other_zeroed_ = .true.
         if (present(other_zeroed)) other_zeroed_ = other_zeroed

         !> Problem size
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Check sizes
         if (n < 1 .or. lda < 1 .or. lda < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a(m,n)=', [lda,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
         ! Compute factorization
         call potrf(triangle,n,a,lda,info)
         call handle_potrf_info(info,triangle,lda,n,err0)
                  
         ! Zero-out the unused part of matrix A
         clean_unused: if (other_zeroed_ .and. err0%ok()) then
             if (lower_) then
                forall (j=2:n) a(:j - 1,j) = zero
             else
                forall (j=1:n - 1) a(j + 1:,j) = zero
             end if
         end if clean_unused

         ! Process output and return
         call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_q_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine stdlib_linalg_q_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(qp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call stdlib_linalg_q_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call linalg_error_handling(err0,err)
         
     end subroutine stdlib_linalg_q_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine stdlib_linalg_c_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,i,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         complex(sp),parameter :: zero = 0.0_sp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         lower_ = .true.
         if (present(lower)) lower_ = lower
         triangle = merge('L','U',lower_)
         
         !> Check if the unused half of the return matrix should be zeroed out (default).
         !> Otherwise it is unused and will contain garbage.
         other_zeroed_ = .true.
         if (present(other_zeroed)) other_zeroed_ = other_zeroed

         !> Problem size
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Check sizes
         if (n < 1 .or. lda < 1 .or. lda < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a(m,n)=', [lda,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
         ! Compute factorization
         call potrf(triangle,n,a,lda,info)
         call handle_potrf_info(info,triangle,lda,n,err0)
                  
         ! Zero-out the unused part of matrix A
         clean_unused: if (other_zeroed_ .and. err0%ok()) then
             if (lower_) then
                forall (j=2:n) a(:j - 1,j) = zero
             else
                forall (j=1:n - 1) a(j + 1:,j) = zero
             end if
         end if clean_unused

         ! Process output and return
         call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_c_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine stdlib_linalg_c_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(sp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call stdlib_linalg_c_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call linalg_error_handling(err0,err)
         
     end subroutine stdlib_linalg_c_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine stdlib_linalg_z_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,i,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         complex(dp),parameter :: zero = 0.0_dp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         lower_ = .true.
         if (present(lower)) lower_ = lower
         triangle = merge('L','U',lower_)
         
         !> Check if the unused half of the return matrix should be zeroed out (default).
         !> Otherwise it is unused and will contain garbage.
         other_zeroed_ = .true.
         if (present(other_zeroed)) other_zeroed_ = other_zeroed

         !> Problem size
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Check sizes
         if (n < 1 .or. lda < 1 .or. lda < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a(m,n)=', [lda,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
         ! Compute factorization
         call potrf(triangle,n,a,lda,info)
         call handle_potrf_info(info,triangle,lda,n,err0)
                  
         ! Zero-out the unused part of matrix A
         clean_unused: if (other_zeroed_ .and. err0%ok()) then
             if (lower_) then
                forall (j=2:n) a(:j - 1,j) = zero
             else
                forall (j=1:n - 1) a(j + 1:,j) = zero
             end if
         end if clean_unused

         ! Process output and return
         call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_z_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine stdlib_linalg_z_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(dp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call stdlib_linalg_z_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call linalg_error_handling(err0,err)
         
     end subroutine stdlib_linalg_z_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine stdlib_linalg_w_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,info,i,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         complex(qp),parameter :: zero = 0.0_qp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         lower_ = .true.
         if (present(lower)) lower_ = lower
         triangle = merge('L','U',lower_)
         
         !> Check if the unused half of the return matrix should be zeroed out (default).
         !> Otherwise it is unused and will contain garbage.
         other_zeroed_ = .true.
         if (present(other_zeroed)) other_zeroed_ = other_zeroed

         !> Problem size
         lda = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)

         ! Check sizes
         if (n < 1 .or. lda < 1 .or. lda < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a(m,n)=', [lda,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
         ! Compute factorization
         call potrf(triangle,n,a,lda,info)
         call handle_potrf_info(info,triangle,lda,n,err0)
                  
         ! Zero-out the unused part of matrix A
         clean_unused: if (other_zeroed_ .and. err0%ok()) then
             if (lower_) then
                forall (j=2:n) a(:j - 1,j) = zero
             else
                forall (j=1:n - 1) a(j + 1:,j) = zero
             end if
         end if clean_unused

         ! Process output and return
         call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_w_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine stdlib_linalg_w_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(qp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix sizes: a=', [lda,n],' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call stdlib_linalg_w_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call linalg_error_handling(err0,err)
         
     end subroutine stdlib_linalg_w_cholesky

end module stdlib_linalg_cholesky
