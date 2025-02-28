!> Cholesky factorization of a matrix, based on LAPACK [POTRF](@ref la_lapack::potrf) functions
module la_cholesky
     use la_constants
     use la_lapack, only: potrf
     use la_state_type, only: la_state, LINALG_ERROR, LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR
     implicit none(type, external)
     private
     
     public :: chol
     public :: cholesky
     
     character(*),parameter :: this = 'cholesky'
     
     !> @brief Computes the Cholesky factorization \f$ A = L \cdot L^T \f$, or \f$ A = U^T \cdot U \f$. 
     !!
     !! ### Summary 
     !! Pure subroutine interface for computing the Cholesky triangular factors. 
     !!
     !! ### Description
     !! 
     !! This interface provides methods for computing the lower- or upper-triangular matrix from the 
     !! Cholesky factorization of a `real` symmetric or `complex` Hermitian matrix.
     !! Supported data types include `real` and `complex`.    
     !! The factorization is computed in-place if only one matrix argument is present; or returned into 
     !! a second matrix argument, if present. The `lower` `logical` flag allows to select between upper or 
     !! lower factorization; the `other_zeroed` optional `logical` flag allows to choose whether the unused
     !! part of the triangular matrix should be filled with zeroes.
     !! 
     !! @note The solution is based on LAPACK's [POTRF](@ref la_lapack::potrf) methods.
     !!
     !! @param[in,out] a The input matrix of size \f$ [m, n] \f$. Overwritten with the Cholesky factor.
     !! @param[out] c (Optional) The output matrix of size \f$ [n, n] \f$, containing the Cholesky factor.
     !! @param[in] lower (Optional) Logical flag indicating whether the lower or upper triangular factor is required.
     !!                  Default: lower triangular factor.
     !! @param[in] other_zeroed (Optional) Logical flag indicating whether the unused half of the matrix should be zeroed.
     !!                         Default: true.
     !! @param[out] err (Optional) State return flag. If not provided, the function will stop on error.
     !!
     interface cholesky
        module procedure la_s_cholesky_inplace
        module procedure la_s_cholesky
        module procedure la_d_cholesky_inplace
        module procedure la_d_cholesky
        module procedure la_q_cholesky_inplace
        module procedure la_q_cholesky
        module procedure la_c_cholesky_inplace
        module procedure la_c_cholesky
        module procedure la_z_cholesky_inplace
        module procedure la_z_cholesky
        module procedure la_w_cholesky_inplace
        module procedure la_w_cholesky
     end interface cholesky         
        
     !> @brief Computes the Cholesky factorization \f$ A = L \cdot L^T \f$, or \f$ A = U^T \cdot U \f$. 
     !!
     !! ### Summary 
     !! Pure function interface for computing the Cholesky triangular factors. 
     !!
     !! ### Description
     !! 
     !! This interface provides methods for computing the lower- or upper-triangular matrix from the 
     !! Cholesky factorization of a `real` symmetric or `complex` Hermitian matrix.
     !! Supported data types include `real` and `complex`.    
     !! The function returns the Cholesky factor as a separate matrix.
     !! The `lower` `logical` flag allows selection between upper or lower factorization,
     !! and the `other_zeroed` optional `logical` flag determines whether the unused
     !! part of the triangular matrix should be zeroed.
     !! 
     !! @note The solution is based on LAPACK's [POTRF](@ref la_lapack::potrf) methods.
     !!
     !! @param[in] a The input matrix of size \f$ [n, n] \f$.
     !! @param[in] lower (Optional) Logical flag indicating whether the lower or upper triangular factor is required.
     !!                  Default: lower triangular factor.
     !! @param[in] other_zeroed (Optional) Logical flag indicating whether the unused half of the matrix should be zeroed.
     !!                         Default: true.
     !! @return c The output matrix of size \f$ [n, n] \f$, containing the Cholesky factor.
     !!
     interface chol
        module procedure la_s_cholesky_fun
        module procedure la_d_cholesky_fun
        module procedure la_q_cholesky_fun
        module procedure la_c_cholesky_fun
        module procedure la_z_cholesky_fun
        module procedure la_w_cholesky_fun
     end interface chol     
     
     

     contains

     elemental subroutine handle_potrf_info(info,triangle,lda,n,err)
         character,intent(in) :: triangle
         integer(ilp),intent(in) :: info,lda,n
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
            case (-4)
                err = la_state(this,LINALG_VALUE_ERROR,'invalid lda=',lda,': is < n = ',n)
            case (1:)
                err = la_state(this,LINALG_ERROR,'cannot complete factorization:',info, &
                                   '-th order leading minor is not positive definite')
            case default
                err = la_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_potrf_info
     
     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine la_s_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         real(sp),parameter :: zero = 0.0_sp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true.
         end if
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
             
             err0 = la_state(this,LINALG_VALUE_ERROR, &
                                      'invalid matrix size: a(m,n)=', [lda,n])
            
         else
         
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
         
         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_s_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine la_s_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(sp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                                     'invalid matrix sizes: a=', [lda,n], &
                                                          ' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call la_s_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call err0%handle(err)
         
     end subroutine la_s_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! Function interface
     pure function la_s_cholesky_fun(a,lower,other_zeroed) result(c)
         !> Input matrix a[n,n]
         real(sp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(sp) :: c(size(a,1),size(a,2))
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         
         c = a
         
         call la_s_cholesky_inplace(c,lower,other_zeroed)

     end function la_s_cholesky_fun

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine la_d_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         real(dp),parameter :: zero = 0.0_dp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true.
         end if
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
             
             err0 = la_state(this,LINALG_VALUE_ERROR, &
                                      'invalid matrix size: a(m,n)=', [lda,n])
            
         else
         
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
         
         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_d_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine la_d_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(dp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                                     'invalid matrix sizes: a=', [lda,n], &
                                                          ' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call la_d_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call err0%handle(err)
         
     end subroutine la_d_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! Function interface
     pure function la_d_cholesky_fun(a,lower,other_zeroed) result(c)
         !> Input matrix a[n,n]
         real(dp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(dp) :: c(size(a,1),size(a,2))
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         
         c = a
         
         call la_d_cholesky_inplace(c,lower,other_zeroed)

     end function la_d_cholesky_fun

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine la_q_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         real(qp),parameter :: zero = 0.0_qp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true.
         end if
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
             
             err0 = la_state(this,LINALG_VALUE_ERROR, &
                                      'invalid matrix size: a(m,n)=', [lda,n])
            
         else
         
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
         
         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_q_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine la_q_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(qp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                                     'invalid matrix sizes: a=', [lda,n], &
                                                          ' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call la_q_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call err0%handle(err)
         
     end subroutine la_q_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! Function interface
     pure function la_q_cholesky_fun(a,lower,other_zeroed) result(c)
         !> Input matrix a[n,n]
         real(qp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         real(qp) :: c(size(a,1),size(a,2))
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         
         c = a
         
         call la_q_cholesky_inplace(c,lower,other_zeroed)

     end function la_q_cholesky_fun

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine la_c_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         complex(sp),parameter :: zero = 0.0_sp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true.
         end if
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
             
             err0 = la_state(this,LINALG_VALUE_ERROR, &
                                      'invalid matrix size: a(m,n)=', [lda,n])
            
         else
         
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
         
         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_c_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine la_c_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(sp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                                     'invalid matrix sizes: a=', [lda,n], &
                                                          ' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call la_c_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call err0%handle(err)
         
     end subroutine la_c_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! Function interface
     pure function la_c_cholesky_fun(a,lower,other_zeroed) result(c)
         !> Input matrix a[n,n]
         complex(sp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(sp) :: c(size(a,1),size(a,2))
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         
         c = a
         
         call la_c_cholesky_inplace(c,lower,other_zeroed)

     end function la_c_cholesky_fun

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine la_z_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         complex(dp),parameter :: zero = 0.0_dp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true.
         end if
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
             
             err0 = la_state(this,LINALG_VALUE_ERROR, &
                                      'invalid matrix size: a(m,n)=', [lda,n])
            
         else
         
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
         
         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_z_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine la_z_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(dp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                                     'invalid matrix sizes: a=', [lda,n], &
                                                          ' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call la_z_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call err0%handle(err)
         
     end subroutine la_z_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! Function interface
     pure function la_z_cholesky_fun(a,lower,other_zeroed) result(c)
         !> Input matrix a[n,n]
         complex(dp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(dp) :: c(size(a,1),size(a,2))
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         
         c = a
         
         call la_z_cholesky_inplace(c,lower,other_zeroed)

     end function la_z_cholesky_fun

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned in-place, overwriting matrix A
     pure subroutine la_w_cholesky_inplace(a,lower,other_zeroed,err)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err

         !> Local variables
         type(la_state) :: err0
         integer(ilp) :: lda,n,info,j
         logical(lk) :: lower_,other_zeroed_
         character :: triangle
         complex(qp),parameter :: zero = 0.0_qp
         
         !> Check if the lower or upper factor is required.
         !> Default: use lower factor
         if (present(lower)) then
            lower_ = lower
         else
            lower_ = .true.
         end if
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
             
             err0 = la_state(this,LINALG_VALUE_ERROR, &
                                      'invalid matrix size: a(m,n)=', [lda,n])
            
         else
         
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
         
         end if

         ! Process output and return
         call err0%handle(err)

     end subroutine la_w_cholesky_inplace

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! The factorization is returned as a separate matrix
     pure subroutine la_w_cholesky(a,c,lower,other_zeroed,err)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(qp),intent(out) :: c(:,:)
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),optional,intent(out) :: err
         
         type(la_state) :: err0
         integer(ilp) :: lda,n,ldc,nc
         
         ! Check C sizes
         lda = size(a,1,kind=ilp)
           n = size(a,2,kind=ilp)
         ldc = size(c,1,kind=ilp)
          nc = size(c,2,kind=ilp)
          
         if (lda < 1 .or. n < 1 .or. lda < n .or. ldc < n .or. nc < n) then
            
            err0 = la_state(this,LINALG_VALUE_ERROR, &
                                     'invalid matrix sizes: a=', [lda,n], &
                                                          ' c=', [ldc,nc])
        
         else
            
            ! Copy data in
            c(:n,:n) = a(:n,:n)
            
            ! Get cholesky factors
            call la_w_cholesky_inplace(c,lower,other_zeroed,err0)
            
         end if
         
         ! Process output and return
         call err0%handle(err)
         
     end subroutine la_w_cholesky

     ! Compute the Cholesky factorization of a symmetric / Hermitian matrix, A = L*L^T = U^T*U.
     ! Function interface
     pure function la_w_cholesky_fun(a,lower,other_zeroed) result(c)
         !> Input matrix a[n,n]
         complex(qp),intent(in) :: a(:,:)
         !> Output matrix with Cholesky factors c[n,n]
         complex(qp) :: c(size(a,1),size(a,2))
         !> [optional] is the lower or upper triangular factor required? Default = lower
         logical(lk),optional,intent(in) :: lower
         !> [optional] should the unused half of the return matrix be zeroed out? Default: yes
         logical(lk),optional,intent(in) :: other_zeroed
         
         c = a
         
         call la_w_cholesky_inplace(c,lower,other_zeroed)

     end function la_w_cholesky_fun

end module la_cholesky
