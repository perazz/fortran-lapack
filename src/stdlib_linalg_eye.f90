
! Return a 2-D matrix with ones on the diagonal and zeros everywhere else
module stdlib_linalg_eye
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Function interface
     public :: eye

     ! Numpy: eye(N, M=None, k=0, dtype=<class 'float'>, order='C', *, device=None, like=None)
     ! Numpy: identity(n, dtype=None, *, like=None) --> square matrices only
     ! Scipy: eye(m, n=None, k=0, dtype=<class 'float'>, format=None) --> sparse only
     ! IMSL:  EYE(N)

     ! Function interface
     interface eye
        module procedure stdlib_linalg_eye_s
        module procedure stdlib_linalg_eye_d
        module procedure stdlib_linalg_eye_q
        module procedure stdlib_linalg_eye_c
        module procedure stdlib_linalg_eye_z
        module procedure stdlib_linalg_eye_w
     end interface eye

     contains

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_s(m,n,dtype,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(sp),intent(in) :: dtype
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Return matrix
         real(sp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(linalg_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            goto 1
         end if

         ! Allocate array
         allocate (eye(m,cols))

         !> Empty matrix
         if (min(m,cols) <= 0) return

         !> Fill data
         do concurrent(i=1:m,j=1:cols)
            eye(i,j) = merge(1.0_sp,0.0_sp,i == j)
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_eye_s

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_d(m,n,dtype,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(dp),optional,intent(in) :: dtype
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Return matrix
         real(dp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(linalg_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            goto 1
         end if

         ! Allocate array
         allocate (eye(m,cols))

         !> Empty matrix
         if (min(m,cols) <= 0) return

         !> Fill data
         do concurrent(i=1:m,j=1:cols)
            eye(i,j) = merge(1.0_dp,0.0_dp,i == j)
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_eye_d

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_q(m,n,dtype,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(qp),intent(in) :: dtype
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Return matrix
         real(qp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(linalg_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            goto 1
         end if

         ! Allocate array
         allocate (eye(m,cols))

         !> Empty matrix
         if (min(m,cols) <= 0) return

         !> Fill data
         do concurrent(i=1:m,j=1:cols)
            eye(i,j) = merge(1.0_qp,0.0_qp,i == j)
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_eye_q

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_c(m,n,dtype,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(sp),intent(in) :: dtype
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Return matrix
         complex(sp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(linalg_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            goto 1
         end if

         ! Allocate array
         allocate (eye(m,cols))

         !> Empty matrix
         if (min(m,cols) <= 0) return

         !> Fill data
         do concurrent(i=1:m,j=1:cols)
            eye(i,j) = merge(1.0_sp,0.0_sp,i == j)
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_eye_c

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_z(m,n,dtype,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(dp),intent(in) :: dtype
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Return matrix
         complex(dp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(linalg_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            goto 1
         end if

         ! Allocate array
         allocate (eye(m,cols))

         !> Empty matrix
         if (min(m,cols) <= 0) return

         !> Fill data
         do concurrent(i=1:m,j=1:cols)
            eye(i,j) = merge(1.0_dp,0.0_dp,i == j)
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_eye_z

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_w(m,n,dtype,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(qp),intent(in) :: dtype
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),optional,intent(out) :: err
         !> Return matrix
         complex(qp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(linalg_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            goto 1
         end if

         ! Allocate array
         allocate (eye(m,cols))

         !> Empty matrix
         if (min(m,cols) <= 0) return

         !> Fill data
         do concurrent(i=1:m,j=1:cols)
            eye(i,j) = merge(1.0_qp,0.0_qp,i == j)
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_eye_w

end module stdlib_linalg_eye
