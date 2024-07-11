! Cholesky factorization of a matrix, based on LAPACK *POTRF functions
module stdlib_linalg_norms
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack,only:lange
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     character(*),parameter :: this = 'norm'
     
     !> List of internal
     integer(ilp),parameter :: NORM_ONE = 1
     integer(ilp),parameter :: NORM_TWO = 2
     integer(ilp),parameter :: NORM_INF = huge(0_ilp) ! infinity norm

     contains
     
     pure subroutine stdlib_linalg_matrix_norm_s(a,order,nrm,err) result(nrm)
         !> Input matrix a[m,n]
         real(sp),intent(in),target :: a(:,:)
         !> Order of the matrix norm being computed.
         integer,intent(in) :: order
         !> Norm of the matrix.
         real(sp),intent(out) :: nrm
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: m,n
         
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         
         ! Initialize norm to zero
         nrm = 0.0_sp
         
         if (.not. (m > 0 .and. n > 0)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
     end function stdlib_linalg_vector_norm_s

     pure subroutine stdlib_linalg_matrix_norm_d(a,order,nrm,err) result(nrm)
         !> Input matrix a[m,n]
         real(dp),intent(in),target :: a(:,:)
         !> Order of the matrix norm being computed.
         integer,intent(in) :: order
         !> Norm of the matrix.
         real(dp),intent(out) :: nrm
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: m,n
         
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         
         ! Initialize norm to zero
         nrm = 0.0_dp
         
         if (.not. (m > 0 .and. n > 0)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
     end function stdlib_linalg_vector_norm_d

     pure subroutine stdlib_linalg_matrix_norm_q(a,order,nrm,err) result(nrm)
         !> Input matrix a[m,n]
         real(qp),intent(in),target :: a(:,:)
         !> Order of the matrix norm being computed.
         integer,intent(in) :: order
         !> Norm of the matrix.
         real(qp),intent(out) :: nrm
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: m,n
         
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         
         ! Initialize norm to zero
         nrm = 0.0_qp
         
         if (.not. (m > 0 .and. n > 0)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
     end function stdlib_linalg_vector_norm_q

     pure subroutine stdlib_linalg_matrix_norm_c(a,order,nrm,err) result(nrm)
         !> Input matrix a[m,n]
         complex(sp),intent(in),target :: a(:,:)
         !> Order of the matrix norm being computed.
         integer,intent(in) :: order
         !> Norm of the matrix.
         real(sp),intent(out) :: nrm
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: m,n
         
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         
         ! Initialize norm to zero
         nrm = 0.0_sp
         
         if (.not. (m > 0 .and. n > 0)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
     end function stdlib_linalg_vector_norm_c

     pure subroutine stdlib_linalg_matrix_norm_z(a,order,nrm,err) result(nrm)
         !> Input matrix a[m,n]
         complex(dp),intent(in),target :: a(:,:)
         !> Order of the matrix norm being computed.
         integer,intent(in) :: order
         !> Norm of the matrix.
         real(dp),intent(out) :: nrm
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: m,n
         
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         
         ! Initialize norm to zero
         nrm = 0.0_dp
         
         if (.not. (m > 0 .and. n > 0)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
     end function stdlib_linalg_vector_norm_z

     pure subroutine stdlib_linalg_matrix_norm_w(a,order,nrm,err) result(nrm)
         !> Input matrix a[m,n]
         complex(qp),intent(in),target :: a(:,:)
         !> Order of the matrix norm being computed.
         integer,intent(in) :: order
         !> Norm of the matrix.
         real(qp),intent(out) :: nrm
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         
         type(linalg_state) :: err0
         integer(ilp) :: m,n
         
         m = size(a,1,kind=ilp)
         n = size(a,2,kind=ilp)
         
         ! Initialize norm to zero
         nrm = 0.0_qp
         
         if (.not. (m > 0 .and. n > 0)) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size: a=', [m,n])
            call linalg_error_handling(err0,err)
            return
         end if
         
     end function stdlib_linalg_vector_norm_w

end module stdlib_linalg_norms
