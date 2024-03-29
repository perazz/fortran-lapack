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
     public :: diag

     ! Numpy: eye(N, M=None, k=0, dtype=<class 'float'>, order='C', *, device=None, like=None)
     ! Numpy: identity(n, dtype=None, *, like=None) --> square matrices only
     ! Scipy: eye(m, n=None, k=0, dtype=<class 'float'>, format=None) --> sparse only
     ! IMSL:  EYE(N)

     ! Identity interface
     interface eye
        module procedure stdlib_linalg_eye_s
        module procedure stdlib_linalg_eye_d
        module procedure stdlib_linalg_eye_q
        module procedure stdlib_linalg_eye_c
        module procedure stdlib_linalg_eye_z
        module procedure stdlib_linalg_eye_w
        module procedure stdlib_linalg_eye_s_errhandle
        module procedure stdlib_linalg_eye_d_errhandle
        module procedure stdlib_linalg_eye_q_errhandle
        module procedure stdlib_linalg_eye_c_errhandle
        module procedure stdlib_linalg_eye_z_errhandle
        module procedure stdlib_linalg_eye_w_errhandle
     end interface eye

     ! Diagonal matrix interface
     interface diag
        module procedure stdlib_linalg_diag_s_from_scalar
        module procedure stdlib_linalg_diag_s_from_array
        module procedure stdlib_linalg_diag_d_from_scalar
        module procedure stdlib_linalg_diag_d_from_array
        module procedure stdlib_linalg_diag_q_from_scalar
        module procedure stdlib_linalg_diag_q_from_array
        module procedure stdlib_linalg_diag_c_from_scalar
        module procedure stdlib_linalg_diag_c_from_array
        module procedure stdlib_linalg_diag_z_from_scalar
        module procedure stdlib_linalg_diag_z_from_array
        module procedure stdlib_linalg_diag_w_from_scalar
        module procedure stdlib_linalg_diag_w_from_array
        module procedure stdlib_linalg_diag_s_errhandle_from_scalar
        module procedure stdlib_linalg_diag_s_errhandle_from_array
        module procedure stdlib_linalg_diag_d_errhandle_from_scalar
        module procedure stdlib_linalg_diag_d_errhandle_from_array
        module procedure stdlib_linalg_diag_q_errhandle_from_scalar
        module procedure stdlib_linalg_diag_q_errhandle_from_array
        module procedure stdlib_linalg_diag_c_errhandle_from_scalar
        module procedure stdlib_linalg_diag_c_errhandle_from_array
        module procedure stdlib_linalg_diag_z_errhandle_from_scalar
        module procedure stdlib_linalg_diag_z_errhandle_from_array
        module procedure stdlib_linalg_diag_w_errhandle_from_scalar
        module procedure stdlib_linalg_diag_w_errhandle_from_array
     end interface diag

     contains

     ! Return diagonal eye matrix of size N
     pure function stdlib_linalg_eye_s(m,n,mold) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(sp),intent(in) :: mold
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
            allocate (eye(0,0))
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
1        call linalg_error_handling(err0)

     end function stdlib_linalg_eye_s

     ! Return diagonal eye matrix of size N
     pure function stdlib_linalg_eye_d(m,n,mold) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(dp),optional,intent(in) :: mold
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
            allocate (eye(0,0))
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
1        call linalg_error_handling(err0)

     end function stdlib_linalg_eye_d

     ! Return diagonal eye matrix of size N
     pure function stdlib_linalg_eye_q(m,n,mold) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(qp),intent(in) :: mold
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
            allocate (eye(0,0))
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
1        call linalg_error_handling(err0)

     end function stdlib_linalg_eye_q

     ! Return diagonal eye matrix of size N
     pure function stdlib_linalg_eye_c(m,n,mold) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(sp),intent(in) :: mold
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
            allocate (eye(0,0))
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
1        call linalg_error_handling(err0)

     end function stdlib_linalg_eye_c

     ! Return diagonal eye matrix of size N
     pure function stdlib_linalg_eye_z(m,n,mold) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(dp),intent(in) :: mold
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
            allocate (eye(0,0))
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
1        call linalg_error_handling(err0)

     end function stdlib_linalg_eye_z

     ! Return diagonal eye matrix of size N
     pure function stdlib_linalg_eye_w(m,n,mold) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(qp),intent(in) :: mold
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
            allocate (eye(0,0))
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
1        call linalg_error_handling(err0)

     end function stdlib_linalg_eye_w

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function stdlib_linalg_diag_s_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(sp),intent(in) :: source
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_s_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function stdlib_linalg_diag_s_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(sp),intent(in) :: source(:)
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_s_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function stdlib_linalg_diag_d_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(dp),intent(in) :: source
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_d_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function stdlib_linalg_diag_d_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(dp),intent(in) :: source(:)
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_d_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function stdlib_linalg_diag_q_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(qp),intent(in) :: source
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_q_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function stdlib_linalg_diag_q_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(qp),intent(in) :: source(:)
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_q_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function stdlib_linalg_diag_c_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(sp),intent(in) :: source
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_c_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function stdlib_linalg_diag_c_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(sp),intent(in) :: source(:)
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_c_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function stdlib_linalg_diag_z_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(dp),intent(in) :: source
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_z_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function stdlib_linalg_diag_z_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(dp),intent(in) :: source(:)
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_z_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function stdlib_linalg_diag_w_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(qp),intent(in) :: source
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_w_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function stdlib_linalg_diag_w_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(qp),intent(in) :: source(:)
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0)

     end function stdlib_linalg_diag_w_from_array

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_s_errhandle(m,n,mold,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(sp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
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
            allocate (eye(0,0))
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

     end function stdlib_linalg_eye_s_errhandle

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_d_errhandle(m,n,mold,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(dp),optional,intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
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
            allocate (eye(0,0))
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

     end function stdlib_linalg_eye_d_errhandle

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_q_errhandle(m,n,mold,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         real(qp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
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
            allocate (eye(0,0))
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

     end function stdlib_linalg_eye_q_errhandle

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_c_errhandle(m,n,mold,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(sp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
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
            allocate (eye(0,0))
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

     end function stdlib_linalg_eye_c_errhandle

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_z_errhandle(m,n,mold,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(dp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
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
            allocate (eye(0,0))
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

     end function stdlib_linalg_eye_z_errhandle

     ! Return diagonal eye matrix of size N
     function stdlib_linalg_eye_w_errhandle(m,n,mold,err) result(eye)
         !> Number of rows
         integer(ilp),intent(in) :: m
         !> Number of columns (optional)
         integer(ilp),optional,intent(in) :: n
         !> Datatype. Used to define the return type. Defaults to real(real64)
         complex(qp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
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
            allocate (eye(0,0))
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

     end function stdlib_linalg_eye_w_errhandle

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function stdlib_linalg_diag_s_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(sp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_s_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function stdlib_linalg_diag_s_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(sp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_s_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function stdlib_linalg_diag_d_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(dp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_d_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function stdlib_linalg_diag_d_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(dp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_d_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function stdlib_linalg_diag_q_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(qp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_q_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function stdlib_linalg_diag_q_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(qp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_q_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function stdlib_linalg_diag_c_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(sp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_c_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function stdlib_linalg_diag_c_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(sp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_sp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_c_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function stdlib_linalg_diag_z_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(dp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_z_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function stdlib_linalg_diag_z_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(dp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_dp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_z_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function stdlib_linalg_diag_w_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(qp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_w_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function stdlib_linalg_diag_w_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(qp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(linalg_state),intent(out) :: err
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(linalg_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            goto 1
         end if

         ! Allocate array
         allocate (diag(n,n))

         !> Empty matrix
         if (n <= 0) return

         !> Fill data
         do concurrent(i=1:n,j=1:n)
            if (i == j) then
               diag(i,j) = source(i)
            else
               diag(i,j) = 0.0_qp
            end if
         end do

         ! Process output and return
1        call linalg_error_handling(err0,err)

     end function stdlib_linalg_diag_w_errhandle_from_array

end module stdlib_linalg_eye
