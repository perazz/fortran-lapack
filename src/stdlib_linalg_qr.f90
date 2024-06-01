module stdlib_linalg_qr
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> QR factorization of a matrix
     public :: qr

     ! Scipy: solve(a, b, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True, assume_a='gen', transposed=False)[source]#
     ! IMSL: lu_solve(a, b, transpose=False)

     interface qr
        module procedure stdlib_linalg_s_qr
        module procedure stdlib_linalg_d_qr
        module procedure stdlib_linalg_q_qr
        module procedure stdlib_linalg_c_qr
        module procedure stdlib_linalg_z_qr
        module procedure stdlib_linalg_w_qr
     end interface qr
     
     character(*),parameter :: this = 'qr'

     contains
     
     elemental subroutine handle_gesv_info(info,m,n,nrhs,err)
         integer(ilp),intent(in) :: info,m,n,nrhs
         type(linalg_state),intent(out) :: err

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=', [m,n])
            case (-7)
                err = linalg_state(this,LINALG_ERROR,'invalid matrix size a=', [m,n])
            case (1:)
                err = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gesv_info

     ! Get workspace size for QR operations
     pure subroutine get_qr_s_workspace(a,lwork,err)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Minimum workspace size for both operations
         integer(ilp),intent(out) :: lwork
         !> State return flag. Returns an error if the query failed
         type(linalg_state),intent(out) :: err
         
         integer(ilp) :: m,n,k,info,lwork_qr,lwork_ord
         real(sp) :: work_dummy(1),tau_dummy(1)
         
         lwork = -1_ilp
         
         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         
         ! QR space
         lwork_qr = -1_ilp
         call geqrf(m,n,a,m,tau_dummy,work_dummy,lwork_qr,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR factorization workspace returned info=',info)
             return
         else
             lwork_qr = ceiling(real(work_dummy(1),kind=sp),kind=ilp)
         end if
         
         ! Ordering space
         lwork_ord = -1_ilp
         call orgqr &
              (m,n,k,a,m,tau_dummy,work_dummy,lwork_ord,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR ordering workspace returned info=',info)
             return
         else
             lwork_ord = ceiling(real(work_dummy(1),kind=sp),kind=ilp)
         end if
         
         ! Pick the largest size, so two operations can be performed with the same allocation
         lwork = max(lwork_qr,lwork_ord)
                  
     end subroutine get_qr_s_workspace
     
     ! Compute the solution to a real system of linear equations A * X = B
     subroutine stdlib_linalg_s_qr(a,q,r,mode,overwrite_a,err)
         !> Input matrix a[m,n]
         real(sp),intent(inout),target :: a(:,:)
         !> Orthogonal matrix Q ([m,m], or [m,k] if reduced)
         real(sp),intent(out),contiguous,target :: q(:,:)
         !> Upper triangular matrix R ([m,n], or [k,n] if reduced)
         real(sp),intent(out),target :: r(:,:)
         !> [optional] Mode: 'reduced' (default), 'complete', 'r', 'raw'
         character(*),optional,intent(in) :: mode
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         character(len=8) :: mode_
         type(linalg_state) :: err0
         integer(ilp) :: m,n,k,q1,q2,r1,r2,lwork,lwork_qr,lwork_ord,info
         logical(lk) :: overwrite_a_
         
         real(sp),pointer :: amat(:,:),tau(:),work(:)

         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         q1 = size(q,1,kind=ilp)
         q2 = size(q,2,kind=ilp)
         r1 = size(r,1,kind=ilp)
         r2 = size(r,2,kind=ilp)

         ! Check sizes
         if (m < 1 .or. n < 1 .or. q1 < m .or. q2 < k .or. r1 < k .or. r2 < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a(m,n)=', [m,n], &
                                                                      ' q(m,m)=', [q1,q2], &
                                                                      ' r(m,n)=', [r1,r2])
            call linalg_error_handling(err0,err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            overwrite_a_ = overwrite_a
         else
            overwrite_a_ = .false._lk
         end if
         
         ! Get mode
         if (present(mode)) then
            mode_ = mode
         else
            mode_ = 'reduced'
         end if

         ! Initialize a matrix temporary
         if (overwrite_a_) then
            amat => a
         else
            allocate (amat(m,n),source=a)
         end if
         
         tau(1:q1*q2) => q

         ! Retrieve workspace size
         call get_qr_s_workspace(a,lwork,err0)

         print *, 'lwork=',lwork
         
         if (err0%ok()) then
                     
             allocate (work(lwork))
             
             ! Compute factorization
             call geqrf(m,n,amat,m,tau,work,lwork,info)
             
             ! Convert K elementary reflectors tau(1:k) -> orthogonal matrix Q
             call orgqr &
                  (m,n,k,amat,m,tau,work,lwork,info)
             
             if (info /= 0) err0 = linalg_state(this,LINALG_VALUE_ERROR,'info=',info)
             
         end if

         if (.not. overwrite_a_) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_s_qr

     ! Get workspace size for QR operations
     pure subroutine get_qr_d_workspace(a,lwork,err)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Minimum workspace size for both operations
         integer(ilp),intent(out) :: lwork
         !> State return flag. Returns an error if the query failed
         type(linalg_state),intent(out) :: err
         
         integer(ilp) :: m,n,k,info,lwork_qr,lwork_ord
         real(dp) :: work_dummy(1),tau_dummy(1)
         
         lwork = -1_ilp
         
         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         
         ! QR space
         lwork_qr = -1_ilp
         call geqrf(m,n,a,m,tau_dummy,work_dummy,lwork_qr,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR factorization workspace returned info=',info)
             return
         else
             lwork_qr = ceiling(real(work_dummy(1),kind=dp),kind=ilp)
         end if
         
         ! Ordering space
         lwork_ord = -1_ilp
         call orgqr &
              (m,n,k,a,m,tau_dummy,work_dummy,lwork_ord,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR ordering workspace returned info=',info)
             return
         else
             lwork_ord = ceiling(real(work_dummy(1),kind=dp),kind=ilp)
         end if
         
         ! Pick the largest size, so two operations can be performed with the same allocation
         lwork = max(lwork_qr,lwork_ord)
                  
     end subroutine get_qr_d_workspace
     
     ! Compute the solution to a real system of linear equations A * X = B
     subroutine stdlib_linalg_d_qr(a,q,r,mode,overwrite_a,err)
         !> Input matrix a[m,n]
         real(dp),intent(inout),target :: a(:,:)
         !> Orthogonal matrix Q ([m,m], or [m,k] if reduced)
         real(dp),intent(out),contiguous,target :: q(:,:)
         !> Upper triangular matrix R ([m,n], or [k,n] if reduced)
         real(dp),intent(out),target :: r(:,:)
         !> [optional] Mode: 'reduced' (default), 'complete', 'r', 'raw'
         character(*),optional,intent(in) :: mode
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         character(len=8) :: mode_
         type(linalg_state) :: err0
         integer(ilp) :: m,n,k,q1,q2,r1,r2,lwork,lwork_qr,lwork_ord,info
         logical(lk) :: overwrite_a_
         
         real(dp),pointer :: amat(:,:),tau(:),work(:)

         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         q1 = size(q,1,kind=ilp)
         q2 = size(q,2,kind=ilp)
         r1 = size(r,1,kind=ilp)
         r2 = size(r,2,kind=ilp)

         ! Check sizes
         if (m < 1 .or. n < 1 .or. q1 < m .or. q2 < k .or. r1 < k .or. r2 < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a(m,n)=', [m,n], &
                                                                      ' q(m,m)=', [q1,q2], &
                                                                      ' r(m,n)=', [r1,r2])
            call linalg_error_handling(err0,err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            overwrite_a_ = overwrite_a
         else
            overwrite_a_ = .false._lk
         end if
         
         ! Get mode
         if (present(mode)) then
            mode_ = mode
         else
            mode_ = 'reduced'
         end if

         ! Initialize a matrix temporary
         if (overwrite_a_) then
            amat => a
         else
            allocate (amat(m,n),source=a)
         end if
         
         tau(1:q1*q2) => q

         ! Retrieve workspace size
         call get_qr_d_workspace(a,lwork,err0)

         print *, 'lwork=',lwork
         
         if (err0%ok()) then
                     
             allocate (work(lwork))
             
             ! Compute factorization
             call geqrf(m,n,amat,m,tau,work,lwork,info)
             
             ! Convert K elementary reflectors tau(1:k) -> orthogonal matrix Q
             call orgqr &
                  (m,n,k,amat,m,tau,work,lwork,info)
             
             if (info /= 0) err0 = linalg_state(this,LINALG_VALUE_ERROR,'info=',info)
             
         end if

         if (.not. overwrite_a_) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_d_qr

     ! Get workspace size for QR operations
     pure subroutine get_qr_q_workspace(a,lwork,err)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Minimum workspace size for both operations
         integer(ilp),intent(out) :: lwork
         !> State return flag. Returns an error if the query failed
         type(linalg_state),intent(out) :: err
         
         integer(ilp) :: m,n,k,info,lwork_qr,lwork_ord
         real(qp) :: work_dummy(1),tau_dummy(1)
         
         lwork = -1_ilp
         
         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         
         ! QR space
         lwork_qr = -1_ilp
         call geqrf(m,n,a,m,tau_dummy,work_dummy,lwork_qr,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR factorization workspace returned info=',info)
             return
         else
             lwork_qr = ceiling(real(work_dummy(1),kind=qp),kind=ilp)
         end if
         
         ! Ordering space
         lwork_ord = -1_ilp
         call orgqr &
              (m,n,k,a,m,tau_dummy,work_dummy,lwork_ord,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR ordering workspace returned info=',info)
             return
         else
             lwork_ord = ceiling(real(work_dummy(1),kind=qp),kind=ilp)
         end if
         
         ! Pick the largest size, so two operations can be performed with the same allocation
         lwork = max(lwork_qr,lwork_ord)
                  
     end subroutine get_qr_q_workspace
     
     ! Compute the solution to a real system of linear equations A * X = B
     subroutine stdlib_linalg_q_qr(a,q,r,mode,overwrite_a,err)
         !> Input matrix a[m,n]
         real(qp),intent(inout),target :: a(:,:)
         !> Orthogonal matrix Q ([m,m], or [m,k] if reduced)
         real(qp),intent(out),contiguous,target :: q(:,:)
         !> Upper triangular matrix R ([m,n], or [k,n] if reduced)
         real(qp),intent(out),target :: r(:,:)
         !> [optional] Mode: 'reduced' (default), 'complete', 'r', 'raw'
         character(*),optional,intent(in) :: mode
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         character(len=8) :: mode_
         type(linalg_state) :: err0
         integer(ilp) :: m,n,k,q1,q2,r1,r2,lwork,lwork_qr,lwork_ord,info
         logical(lk) :: overwrite_a_
         
         real(qp),pointer :: amat(:,:),tau(:),work(:)

         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         q1 = size(q,1,kind=ilp)
         q2 = size(q,2,kind=ilp)
         r1 = size(r,1,kind=ilp)
         r2 = size(r,2,kind=ilp)

         ! Check sizes
         if (m < 1 .or. n < 1 .or. q1 < m .or. q2 < k .or. r1 < k .or. r2 < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a(m,n)=', [m,n], &
                                                                      ' q(m,m)=', [q1,q2], &
                                                                      ' r(m,n)=', [r1,r2])
            call linalg_error_handling(err0,err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            overwrite_a_ = overwrite_a
         else
            overwrite_a_ = .false._lk
         end if
         
         ! Get mode
         if (present(mode)) then
            mode_ = mode
         else
            mode_ = 'reduced'
         end if

         ! Initialize a matrix temporary
         if (overwrite_a_) then
            amat => a
         else
            allocate (amat(m,n),source=a)
         end if
         
         tau(1:q1*q2) => q

         ! Retrieve workspace size
         call get_qr_q_workspace(a,lwork,err0)

         print *, 'lwork=',lwork
         
         if (err0%ok()) then
                     
             allocate (work(lwork))
             
             ! Compute factorization
             call geqrf(m,n,amat,m,tau,work,lwork,info)
             
             ! Convert K elementary reflectors tau(1:k) -> orthogonal matrix Q
             call orgqr &
                  (m,n,k,amat,m,tau,work,lwork,info)
             
             if (info /= 0) err0 = linalg_state(this,LINALG_VALUE_ERROR,'info=',info)
             
         end if

         if (.not. overwrite_a_) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_q_qr

     ! Get workspace size for QR operations
     pure subroutine get_qr_c_workspace(a,lwork,err)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Minimum workspace size for both operations
         integer(ilp),intent(out) :: lwork
         !> State return flag. Returns an error if the query failed
         type(linalg_state),intent(out) :: err
         
         integer(ilp) :: m,n,k,info,lwork_qr,lwork_ord
         complex(sp) :: work_dummy(1),tau_dummy(1)
         
         lwork = -1_ilp
         
         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         
         ! QR space
         lwork_qr = -1_ilp
         call geqrf(m,n,a,m,tau_dummy,work_dummy,lwork_qr,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR factorization workspace returned info=',info)
             return
         else
             lwork_qr = ceiling(real(work_dummy(1),kind=sp),kind=ilp)
         end if
         
         ! Ordering space
         lwork_ord = -1_ilp
         call ungqr &
              (m,n,k,a,m,tau_dummy,work_dummy,lwork_ord,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR ordering workspace returned info=',info)
             return
         else
             lwork_ord = ceiling(real(work_dummy(1),kind=sp),kind=ilp)
         end if
         
         ! Pick the largest size, so two operations can be performed with the same allocation
         lwork = max(lwork_qr,lwork_ord)
                  
     end subroutine get_qr_c_workspace
     
     ! Compute the solution to a real system of linear equations A * X = B
     subroutine stdlib_linalg_c_qr(a,q,r,mode,overwrite_a,err)
         !> Input matrix a[m,n]
         complex(sp),intent(inout),target :: a(:,:)
         !> Orthogonal matrix Q ([m,m], or [m,k] if reduced)
         complex(sp),intent(out),contiguous,target :: q(:,:)
         !> Upper triangular matrix R ([m,n], or [k,n] if reduced)
         complex(sp),intent(out),target :: r(:,:)
         !> [optional] Mode: 'reduced' (default), 'complete', 'r', 'raw'
         character(*),optional,intent(in) :: mode
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         character(len=8) :: mode_
         type(linalg_state) :: err0
         integer(ilp) :: m,n,k,q1,q2,r1,r2,lwork,lwork_qr,lwork_ord,info
         logical(lk) :: overwrite_a_
         
         complex(sp),pointer :: amat(:,:),tau(:),work(:)

         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         q1 = size(q,1,kind=ilp)
         q2 = size(q,2,kind=ilp)
         r1 = size(r,1,kind=ilp)
         r2 = size(r,2,kind=ilp)

         ! Check sizes
         if (m < 1 .or. n < 1 .or. q1 < m .or. q2 < k .or. r1 < k .or. r2 < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a(m,n)=', [m,n], &
                                                                      ' q(m,m)=', [q1,q2], &
                                                                      ' r(m,n)=', [r1,r2])
            call linalg_error_handling(err0,err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            overwrite_a_ = overwrite_a
         else
            overwrite_a_ = .false._lk
         end if
         
         ! Get mode
         if (present(mode)) then
            mode_ = mode
         else
            mode_ = 'reduced'
         end if

         ! Initialize a matrix temporary
         if (overwrite_a_) then
            amat => a
         else
            allocate (amat(m,n),source=a)
         end if
         
         tau(1:q1*q2) => q

         ! Retrieve workspace size
         call get_qr_c_workspace(a,lwork,err0)

         print *, 'lwork=',lwork
         
         if (err0%ok()) then
                     
             allocate (work(lwork))
             
             ! Compute factorization
             call geqrf(m,n,amat,m,tau,work,lwork,info)
             
             ! Convert K elementary reflectors tau(1:k) -> orthogonal matrix Q
             call ungqr &
                  (m,n,k,amat,m,tau,work,lwork,info)
             
             if (info /= 0) err0 = linalg_state(this,LINALG_VALUE_ERROR,'info=',info)
             
         end if

         if (.not. overwrite_a_) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_c_qr

     ! Get workspace size for QR operations
     pure subroutine get_qr_z_workspace(a,lwork,err)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Minimum workspace size for both operations
         integer(ilp),intent(out) :: lwork
         !> State return flag. Returns an error if the query failed
         type(linalg_state),intent(out) :: err
         
         integer(ilp) :: m,n,k,info,lwork_qr,lwork_ord
         complex(dp) :: work_dummy(1),tau_dummy(1)
         
         lwork = -1_ilp
         
         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         
         ! QR space
         lwork_qr = -1_ilp
         call geqrf(m,n,a,m,tau_dummy,work_dummy,lwork_qr,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR factorization workspace returned info=',info)
             return
         else
             lwork_qr = ceiling(real(work_dummy(1),kind=dp),kind=ilp)
         end if
         
         ! Ordering space
         lwork_ord = -1_ilp
         call ungqr &
              (m,n,k,a,m,tau_dummy,work_dummy,lwork_ord,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR ordering workspace returned info=',info)
             return
         else
             lwork_ord = ceiling(real(work_dummy(1),kind=dp),kind=ilp)
         end if
         
         ! Pick the largest size, so two operations can be performed with the same allocation
         lwork = max(lwork_qr,lwork_ord)
                  
     end subroutine get_qr_z_workspace
     
     ! Compute the solution to a real system of linear equations A * X = B
     subroutine stdlib_linalg_z_qr(a,q,r,mode,overwrite_a,err)
         !> Input matrix a[m,n]
         complex(dp),intent(inout),target :: a(:,:)
         !> Orthogonal matrix Q ([m,m], or [m,k] if reduced)
         complex(dp),intent(out),contiguous,target :: q(:,:)
         !> Upper triangular matrix R ([m,n], or [k,n] if reduced)
         complex(dp),intent(out),target :: r(:,:)
         !> [optional] Mode: 'reduced' (default), 'complete', 'r', 'raw'
         character(*),optional,intent(in) :: mode
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         character(len=8) :: mode_
         type(linalg_state) :: err0
         integer(ilp) :: m,n,k,q1,q2,r1,r2,lwork,lwork_qr,lwork_ord,info
         logical(lk) :: overwrite_a_
         
         complex(dp),pointer :: amat(:,:),tau(:),work(:)

         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         q1 = size(q,1,kind=ilp)
         q2 = size(q,2,kind=ilp)
         r1 = size(r,1,kind=ilp)
         r2 = size(r,2,kind=ilp)

         ! Check sizes
         if (m < 1 .or. n < 1 .or. q1 < m .or. q2 < k .or. r1 < k .or. r2 < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a(m,n)=', [m,n], &
                                                                      ' q(m,m)=', [q1,q2], &
                                                                      ' r(m,n)=', [r1,r2])
            call linalg_error_handling(err0,err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            overwrite_a_ = overwrite_a
         else
            overwrite_a_ = .false._lk
         end if
         
         ! Get mode
         if (present(mode)) then
            mode_ = mode
         else
            mode_ = 'reduced'
         end if

         ! Initialize a matrix temporary
         if (overwrite_a_) then
            amat => a
         else
            allocate (amat(m,n),source=a)
         end if
         
         tau(1:q1*q2) => q

         ! Retrieve workspace size
         call get_qr_z_workspace(a,lwork,err0)

         print *, 'lwork=',lwork
         
         if (err0%ok()) then
                     
             allocate (work(lwork))
             
             ! Compute factorization
             call geqrf(m,n,amat,m,tau,work,lwork,info)
             
             ! Convert K elementary reflectors tau(1:k) -> orthogonal matrix Q
             call ungqr &
                  (m,n,k,amat,m,tau,work,lwork,info)
             
             if (info /= 0) err0 = linalg_state(this,LINALG_VALUE_ERROR,'info=',info)
             
         end if

         if (.not. overwrite_a_) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_z_qr

     ! Get workspace size for QR operations
     pure subroutine get_qr_w_workspace(a,lwork,err)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Minimum workspace size for both operations
         integer(ilp),intent(out) :: lwork
         !> State return flag. Returns an error if the query failed
         type(linalg_state),intent(out) :: err
         
         integer(ilp) :: m,n,k,info,lwork_qr,lwork_ord
         complex(qp) :: work_dummy(1),tau_dummy(1)
         
         lwork = -1_ilp
         
         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         
         ! QR space
         lwork_qr = -1_ilp
         call geqrf(m,n,a,m,tau_dummy,work_dummy,lwork_qr,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR factorization workspace returned info=',info)
             return
         else
             lwork_qr = ceiling(real(work_dummy(1),kind=qp),kind=ilp)
         end if
         
         ! Ordering space
         lwork_ord = -1_ilp
         call ungqr &
              (m,n,k,a,m,tau_dummy,work_dummy,lwork_ord,info)
         if (info /= 0) then
             err = linalg_state(this,LINALG_INTERNAL_ERROR,'QR ordering workspace returned info=',info)
             return
         else
             lwork_ord = ceiling(real(work_dummy(1),kind=qp),kind=ilp)
         end if
         
         ! Pick the largest size, so two operations can be performed with the same allocation
         lwork = max(lwork_qr,lwork_ord)
                  
     end subroutine get_qr_w_workspace
     
     ! Compute the solution to a real system of linear equations A * X = B
     subroutine stdlib_linalg_w_qr(a,q,r,mode,overwrite_a,err)
         !> Input matrix a[m,n]
         complex(qp),intent(inout),target :: a(:,:)
         !> Orthogonal matrix Q ([m,m], or [m,k] if reduced)
         complex(qp),intent(out),contiguous,target :: q(:,:)
         !> Upper triangular matrix R ([m,n], or [k,n] if reduced)
         complex(qp),intent(out),target :: r(:,:)
         !> [optional] Mode: 'reduced' (default), 'complete', 'r', 'raw'
         character(*),optional,intent(in) :: mode
         !> [optional] Can A data be overwritten and destroyed?
         logical(lk),optional,intent(in) :: overwrite_a
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err

         !> Local variables
         character(len=8) :: mode_
         type(linalg_state) :: err0
         integer(ilp) :: m,n,k,q1,q2,r1,r2,lwork,lwork_qr,lwork_ord,info
         logical(lk) :: overwrite_a_
         
         complex(qp),pointer :: amat(:,:),tau(:),work(:)

         !> Problem sizes
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         k = min(m,n)
         q1 = size(q,1,kind=ilp)
         q2 = size(q,2,kind=ilp)
         r1 = size(r,1,kind=ilp)
         r2 = size(r,2,kind=ilp)

         ! Check sizes
         if (m < 1 .or. n < 1 .or. q1 < m .or. q2 < k .or. r1 < k .or. r2 < n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a(m,n)=', [m,n], &
                                                                      ' q(m,m)=', [q1,q2], &
                                                                      ' r(m,n)=', [r1,r2])
            call linalg_error_handling(err0,err)
            return
         end if

         ! Can A be overwritten? By default, do not overwrite
         if (present(overwrite_a)) then
            overwrite_a_ = overwrite_a
         else
            overwrite_a_ = .false._lk
         end if
         
         ! Get mode
         if (present(mode)) then
            mode_ = mode
         else
            mode_ = 'reduced'
         end if

         ! Initialize a matrix temporary
         if (overwrite_a_) then
            amat => a
         else
            allocate (amat(m,n),source=a)
         end if
         
         tau(1:q1*q2) => q

         ! Retrieve workspace size
         call get_qr_w_workspace(a,lwork,err0)

         print *, 'lwork=',lwork
         
         if (err0%ok()) then
                     
             allocate (work(lwork))
             
             ! Compute factorization
             call geqrf(m,n,amat,m,tau,work,lwork,info)
             
             ! Convert K elementary reflectors tau(1:k) -> orthogonal matrix Q
             call ungqr &
                  (m,n,k,amat,m,tau,work,lwork,info)
             
             if (info /= 0) err0 = linalg_state(this,LINALG_VALUE_ERROR,'info=',info)
             
         end if

         if (.not. overwrite_a_) deallocate (amat)

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end subroutine stdlib_linalg_w_qr

end module stdlib_linalg_qr
