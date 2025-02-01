!> Identity and diagonal matrices.
module la_eye
     use la_constants
     use la_blas
     use la_lapack
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

    !> @brief Construct an identity matrix of size \f$m \times n\f$.
    !!
    !! This function returns a diagonal identity matrix of size \f$m \times n\f$, where all diagonal elements
    !! are set to 1 and all off-diagonal elements are set to 0. The number of rows and columns can be specified.
    !! If only one parameter is provided, a square matrix of size \f$m \times m\f$ is returned.
    !!
    !! @param[in] m The number of rows of the identity matrix.
    !! @param[in] n (Optional) The number of columns of the identity matrix. If omitted, the matrix is square (\f$m \times m\f$).
    !! @param[in] mold (Optional) Data type to define the return type. Defaults to `real(real64)`.
    !! 
    !! @return The identity matrix with size \f$m \times n\f$.
    !!
    !! @note If the `mold` parameter is omitted, the default type is `real(real64)`. If specified, the return type
    !!       matches the given type.
    !!
    !! @warning Ensure that the matrix dimensions are valid and consistent with the type definition.
    public :: eye

    !> @brief Return a square diagonal matrix with diagonal values.
    !!
    !! This function generates a square diagonal matrix where the diagonal elements are either
    !! equal to the specified scalar value or populated by the input array.
    !! The size of the matrix is determined by the input parameter \f$n\f$ or the size of the input array.
    !!
    !! @param[in] n The size of the square matrix (only used if a scalar is provided for the diagonal).
    !! @param[in] source If a scalar, this value is used to populate the diagonal of the matrix.
    !!                   If an array, the elements of the array are used for the diagonal.
    !! @param[out] err (Optional) State return flag. If not provided, the function will stop on error.
    !! @return The diagonal matrix with size \f$n \times n\f$, where the diagonal elements are populated by \f$source\f$.
    !!
    !! @note If a scalar value is passed, the diagonal elements of the matrix will all be equal to \f$source\f$.
    !!       If an array is passed, its length determines the size of the matrix, and the array elements are placed
    !!       along the diagonal. The `err` parameter is optional. If not requested, the code will stop on error.
    !!       Otherwise, it returns the error state of the function.
    !!
    public :: diag

     interface eye
        module procedure la_eye_d
        module procedure la_eye_q
        module procedure la_eye_c
        module procedure la_eye_z
        module procedure la_eye_w
        module procedure la_eye_s_errhandle
        module procedure la_eye_d_errhandle
        module procedure la_eye_q_errhandle
        module procedure la_eye_c_errhandle
        module procedure la_eye_z_errhandle
        module procedure la_eye_w_errhandle
     end interface eye

     ! Diagonal matrix interface
     interface diag
        module procedure la_diag_s_from_scalar
        module procedure la_diag_s_from_array
        module procedure la_diag_d_from_scalar
        module procedure la_diag_d_from_array
        module procedure la_diag_q_from_scalar
        module procedure la_diag_q_from_array
        module procedure la_diag_c_from_scalar
        module procedure la_diag_c_from_array
        module procedure la_diag_z_from_scalar
        module procedure la_diag_z_from_array
        module procedure la_diag_w_from_scalar
        module procedure la_diag_w_from_array
        module procedure la_diag_s_errhandle_from_scalar
        module procedure la_diag_s_errhandle_from_array
        module procedure la_diag_d_errhandle_from_scalar
        module procedure la_diag_d_errhandle_from_array
        module procedure la_diag_q_errhandle_from_scalar
        module procedure la_diag_q_errhandle_from_array
        module procedure la_diag_c_errhandle_from_scalar
        module procedure la_diag_c_errhandle_from_array
        module procedure la_diag_z_errhandle_from_scalar
        module procedure la_diag_z_errhandle_from_array
        module procedure la_diag_w_errhandle_from_scalar
        module procedure la_diag_w_errhandle_from_array
     end interface diag

     contains

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     pure function la_eye_s(m,n,mold) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         real(sp),intent(in) :: mold
         !> Return matrix
         real(sp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_eye_s

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     pure function la_eye_d(m,n,mold) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         real(dp),optional,intent(in) :: mold
         !> Return matrix
         real(dp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_eye_d

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     pure function la_eye_q(m,n,mold) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         real(qp),intent(in) :: mold
         !> Return matrix
         real(qp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_eye_q

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     pure function la_eye_c(m,n,mold) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         complex(sp),intent(in) :: mold
         !> Return matrix
         complex(sp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_eye_c

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     pure function la_eye_z(m,n,mold) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         complex(dp),intent(in) :: mold
         !> Return matrix
         complex(dp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_eye_z

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     pure function la_eye_w(m,n,mold) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         complex(qp),intent(in) :: mold
         !> Return matrix
         complex(qp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_eye_w

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function la_diag_s_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(sp),intent(in) :: source
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_diag_s_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function la_diag_s_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(sp),intent(in) :: source(:)
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
         call err0%handle()

     end function la_diag_s_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function la_diag_d_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(dp),intent(in) :: source
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_diag_d_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function la_diag_d_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(dp),intent(in) :: source(:)
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
         call err0%handle()

     end function la_diag_d_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function la_diag_q_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(qp),intent(in) :: source
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_diag_q_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function la_diag_q_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(qp),intent(in) :: source(:)
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
         call err0%handle()

     end function la_diag_q_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function la_diag_c_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(sp),intent(in) :: source
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_diag_c_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function la_diag_c_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(sp),intent(in) :: source(:)
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
         call err0%handle()

     end function la_diag_c_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function la_diag_z_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(dp),intent(in) :: source
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_diag_z_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function la_diag_z_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(dp),intent(in) :: source(:)
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
         call err0%handle()

     end function la_diag_z_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     pure function la_diag_w_from_scalar(n,source) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(qp),intent(in) :: source
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
1        call err0%handle()

     end function la_diag_w_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     pure function la_diag_w_from_array(source) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(qp),intent(in) :: source(:)
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle()
            return
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
         call err0%handle()

     end function la_diag_w_from_array

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     function la_eye_s_errhandle(m,n,mold,err) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         real(sp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(sp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_eye_s_errhandle

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     function la_eye_d_errhandle(m,n,mold,err) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         real(dp),optional,intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(dp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_eye_d_errhandle

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     function la_eye_q_errhandle(m,n,mold,err) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         real(qp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(qp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_eye_q_errhandle

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     function la_eye_c_errhandle(m,n,mold,err) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         complex(sp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(sp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_eye_c_errhandle

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     function la_eye_z_errhandle(m,n,mold,err) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         complex(dp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(dp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_eye_z_errhandle

     !> Function to construct an identity matrix of size `m x n`.
     !! This function returns a diagonal identity matrix, with the diagonal elements
     !! equal to 1 and all other elements set to 0.
     function la_eye_w_errhandle(m,n,mold,err) result(eye)
         !> Number of rows of the identity matrix.
         integer(ilp),intent(in) :: m
         !> Number of columns of the identity matrix (optional).
         integer(ilp),optional,intent(in) :: n
         !> Data type, used to define the return type. Defaults to `real(real64)`.
         complex(qp),intent(in) :: mold
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(qp),allocatable :: eye(:,:)

         !> Local variables
         integer(ilp) :: i,j,cols
         type(la_state) :: err0
         character(*),parameter :: this = 'eye'

         !> Determine number of columns
         if (present(n)) then
            cols = n
         else
            cols = m
         end if

         !> Check size
         if (.not. min(m,cols) >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid eye size: eye[',m,',',n,']')
            allocate (eye(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_eye_w_errhandle

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function la_diag_s_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(sp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_diag_s_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function la_diag_s_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(sp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
         call err0%handle(err)

     end function la_diag_s_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function la_diag_d_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(dp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_diag_d_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function la_diag_d_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(dp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
         call err0%handle(err)

     end function la_diag_d_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function la_diag_q_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         real(qp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_diag_q_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function la_diag_q_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(qp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         real(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
         call err0%handle(err)

     end function la_diag_q_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function la_diag_c_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(sp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_diag_c_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function la_diag_c_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(sp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(sp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
         call err0%handle(err)

     end function la_diag_c_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function la_diag_z_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(dp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_diag_z_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function la_diag_z_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(dp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(dp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
         call err0%handle(err)

     end function la_diag_z_errhandle_from_array

     ! Return square diagonal matrix with diagonal values equal to the input scalar
     function la_diag_w_errhandle_from_scalar(n,source,err) result(diag)
         !> Matrix size
         integer(ilp),intent(in) :: n
         !> Scalar diagonal value. Used to define the return type.
         complex(qp),intent(in) :: source
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::scalar'

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid diagonal size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
1        call err0%handle(err)

     end function la_diag_w_errhandle_from_scalar

     ! Construct square diagonal matrix from an array of diagonal values
     function la_diag_w_errhandle_from_array(source,err) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(qp),intent(in) :: source(:)
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err
         !> Return matrix
         complex(qp),allocatable :: diag(:,:)

         !> Local variables
         integer(ilp) :: i,j,n
         type(la_state) :: err0
         character(*),parameter :: this = 'diag::array'

         n = size(source,kind=ilp)

         !> Check size
         if (.not. n >= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid input array size: diag[',n,',',n,']')
            allocate (diag(0,0))
            call err0%handle(err)
            return
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
         call err0%handle(err)

     end function la_diag_w_errhandle_from_array

end module la_eye
