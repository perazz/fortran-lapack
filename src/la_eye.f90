!> Identity and diagonal matrices, matrix trace and elementary matrix products.
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

    !> @brief Return a square diagonal matrix with diagonal values, or extract a diagonal.
    !!
    !! This function generates a square diagonal matrix where the diagonal elements are either
    !! equal to the specified scalar value or populated by the input array.
    !! The size of the matrix is determined by the input parameter \f$n\f$ or the size of the input array.
    !! Given a matrix instead, the function extracts one of its diagonals as a vector.
    !!
    !! @param[in] n The size of the square matrix (only used if a scalar is provided for the diagonal).
    !! @param[in] source If a scalar, this value is used to populate the diagonal of the matrix.
    !!                   If an array, the elements of the array are used for the diagonal.
    !!                   If a matrix, its \f$k\f$-th diagonal is returned as a vector.
    !! @param[in] k (Optional) Index of the diagonal: 0 is the main diagonal, \f$k>0\f$ the \f$k\f$-th
    !!              superdiagonal, \f$k<0\f$ the \f$k\f$-th subdiagonal.
    !! @param[out] err (Optional) State return flag. If not provided, the function will stop on error.
    !! @return The diagonal matrix with size \f$n \times n\f$, where the diagonal elements are populated by \f$source\f$,
    !!         or the requested diagonal of the input matrix.
    !!
    !! @note If a scalar value is passed, the diagonal elements of the matrix will all be equal to \f$source\f$.
    !!       If an array is passed, its length determines the size of the matrix, and the array elements are placed
    !!       along the diagonal. The `err` parameter is optional. If not requested, the code will stop on error.
    !!       Otherwise, it returns the error state of the function.
    !!
    public :: diag

    !> @brief Return the trace of a matrix.
    !!
    !! This function sums the elements of the main diagonal of a matrix. The matrix need not be square:
    !! for a \f$m \times n\f$ matrix, the first \f$\min(m,n)\f$ diagonal elements are summed.
    !!
    !! @param[in] A The input matrix of size \f$ [m,n] \f$.
    !! @return The sum \f$ \sum_i A_{ii} \f$, of the same type and kind as `A`.
    !!
    public :: trace

    !> @brief Return the outer product of two vectors.
    !!
    !! This function computes \f$ u \otimes v \f$, the rank-2 array whose \f$(i,j)\f$ element is
    !! \f$ u_i v_j \f$.
    !!
    !! @param[in] u The first input vector, of size \f$m\f$.
    !! @param[in] v The second input vector, of size \f$n\f$.
    !! @return The \f$ m \times n \f$ outer product matrix.
    !!
    public :: outer_product

    !> @brief Return the cross product of two 3-dimensional vectors.
    !!
    !! This function computes \f$ a \times b \f$, the vector orthogonal to both input vectors.
    !!
    !! @param[in] a The first input vector, of size 3.
    !! @param[in] b The second input vector, of size 3.
    !! @return The size-3 cross product vector.
    !!
    public :: cross_product

    !> @brief Return the Kronecker product of two matrices.
    !!
    !! This function computes \f$ A \otimes B \f$: given \f$ A \f$ of size \f$ m_1 \times n_1 \f$ and
    !! \f$ B \f$ of size \f$ m_2 \times n_2 \f$, the result is the \f$ (m_1 m_2) \times (n_1 n_2) \f$
    !! block matrix whose \f$(i,j)\f$ block is \f$ A_{ij} B \f$.
    !!
    !! @param[in] A The first input matrix, of size \f$ [m_1,n_1] \f$.
    !! @param[in] B The second input matrix, of size \f$ [m_2,n_2] \f$.
    !! @return The \f$ (m_1 m_2) \times (n_1 n_2) \f$ Kronecker product matrix.
    !!
    public :: kronecker_product

    !> @brief Return the Hermitian transpose of a matrix.
    !!
    !! This function returns `conjg(transpose(a))` for a complex matrix and `transpose(a)` for a real one.
    !!
    !! @param[in] a The input matrix of size \f$ [m,n] \f$.
    !! @return The \f$ n \times m \f$ matrix \f$ a^H \f$.
    !!
    public :: hermitian

     ! Identity matrix interface
     interface eye
        module procedure la_eye_s
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
        module procedure la_diag_s_from_array_k
        module procedure la_diag_s_from_matrix
        module procedure la_diag_s_from_matrix_k
        module procedure la_diag_d_from_array_k
        module procedure la_diag_d_from_matrix
        module procedure la_diag_d_from_matrix_k
        module procedure la_diag_q_from_array_k
        module procedure la_diag_q_from_matrix
        module procedure la_diag_q_from_matrix_k
        module procedure la_diag_c_from_array_k
        module procedure la_diag_c_from_matrix
        module procedure la_diag_c_from_matrix_k
        module procedure la_diag_z_from_array_k
        module procedure la_diag_z_from_matrix
        module procedure la_diag_z_from_matrix_k
        module procedure la_diag_w_from_array_k
        module procedure la_diag_w_from_matrix
        module procedure la_diag_w_from_matrix_k
     end interface diag

     ! Matrix trace interface
     interface trace
        module procedure la_trace_s
        module procedure la_trace_d
        module procedure la_trace_q
        module procedure la_trace_c
        module procedure la_trace_z
        module procedure la_trace_w
     end interface trace

     ! Outer product interface
     interface outer_product
        module procedure la_outer_product_s
        module procedure la_outer_product_d
        module procedure la_outer_product_q
        module procedure la_outer_product_c
        module procedure la_outer_product_z
        module procedure la_outer_product_w
     end interface outer_product

     ! Cross product interface
     interface cross_product
        module procedure la_cross_product_s
        module procedure la_cross_product_d
        module procedure la_cross_product_q
        module procedure la_cross_product_c
        module procedure la_cross_product_z
        module procedure la_cross_product_w
     end interface cross_product

     ! Kronecker product interface
     interface kronecker_product
        module procedure la_kronecker_product_s
        module procedure la_kronecker_product_d
        module procedure la_kronecker_product_q
        module procedure la_kronecker_product_c
        module procedure la_kronecker_product_z
        module procedure la_kronecker_product_w
     end interface kronecker_product

     ! Hermitian transpose interface
     interface hermitian
        module procedure la_hermitian_s
        module procedure la_hermitian_d
        module procedure la_hermitian_q
        module procedure la_hermitian_c
        module procedure la_hermitian_z
        module procedure la_hermitian_w
     end interface hermitian

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

     ! Construct a square matrix whose k-th diagonal is filled with the input array
     pure function la_diag_s_from_array_k(source,k) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(sp),intent(in) :: source(:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return matrix
         real(sp) :: diag(size(source,kind=ilp) + abs(k),size(source,kind=ilp) + abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = size(source,kind=ilp)

         diag = 0.0_sp

         if (k > 0) then
            do i = 1,n
               diag(i,k + i) = source(i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i + abs(k),i) = source(i)
            end do
         else
            do i = 1,n
               diag(i,i) = source(i)
            end do
         end if

     end function la_diag_s_from_array_k

     ! Extract the main diagonal of a matrix
     pure function la_diag_s_from_matrix(A) result(diag)
         !> Input matrix. Used to define the return type and size.
         real(sp),intent(in) :: A(:,:)
         !> Return array
         real(sp) :: diag(minval(shape(A,kind=ilp)))

         !> Local variables
         integer(ilp) :: i

         do i = 1,minval(shape(A,kind=ilp))
            diag(i) = A(i,i)
         end do

     end function la_diag_s_from_matrix

     ! Extract the k-th diagonal of a matrix
     pure function la_diag_s_from_matrix_k(A,k) result(diag)
         !> Input matrix. Used to define the return type and size.
         real(sp),intent(in) :: A(:,:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return array
         real(sp) :: diag(minval(shape(A,kind=ilp)) - abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = minval(shape(A,kind=ilp)) - abs(k)

         if (k > 0) then
            do i = 1,n
               diag(i) = A(i,k + i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i) = A(i + abs(k),i)
            end do
         else
            do i = 1,n
               diag(i) = A(i,i)
            end do
         end if

     end function la_diag_s_from_matrix_k

     ! Sum of the main diagonal elements of a matrix
     pure function la_trace_s(A) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)
         !> Trace of the input matrix
         real(sp) :: res

         !> Local variables
         integer(ilp) :: i

         res = 0.0_sp
         do i = 1,minval(shape(A,kind=ilp))
            res = res + A(i,i)
         end do

     end function la_trace_s

     ! Outer product of two vectors
     pure function la_outer_product_s(u,v) result(res)
         !> First input vector
         real(sp),intent(in) :: u(:)
         !> Second input vector
         real(sp),intent(in) :: v(:)
         !> Outer product matrix
         real(sp) :: res(size(u,kind=ilp),size(v,kind=ilp))

         !> Local variables
         integer(ilp) :: col

         do col = 1,size(v,kind=ilp)
            res(:,col) = v(col)*u
         end do

     end function la_outer_product_s

     ! Cross product of two 3-dimensional vectors
     pure function la_cross_product_s(a,b) result(res)
         !> First input vector
         real(sp),intent(in) :: a(3)
         !> Second input vector
         real(sp),intent(in) :: b(3)
         !> Cross product vector
         real(sp) :: res(3)

         res(1) = a(2)*b(3) - a(3)*b(2)
         res(2) = a(3)*b(1) - a(1)*b(3)
         res(3) = a(1)*b(2) - a(2)*b(1)

     end function la_cross_product_s

     ! Kronecker product of two matrices
     pure function la_kronecker_product_s(A,B) result(C)
         !> First input matrix
         real(sp),intent(in) :: A(:,:)
         !> Second input matrix
         real(sp),intent(in) :: B(:,:)
         !> Kronecker product matrix
         real(sp) :: C(size(A,1,kind=ilp)*size(B,1,kind=ilp),size(A,2,kind=ilp)*size(B,2,kind=ilp))

         !> Local variables
         integer(ilp) :: m1,n1,maxM1,maxN1,maxM2,maxN2

         maxM1 = size(A,1,kind=ilp)
         maxN1 = size(A,2,kind=ilp)
         maxM2 = size(B,1,kind=ilp)
         maxN2 = size(B,2,kind=ilp)

         do n1 = 1,maxN1
            do m1 = 1,maxM1
               C((m1 - 1)*maxM2 + 1:m1*maxM2, (n1 - 1)*maxN2 + 1:n1*maxN2) = A(m1,n1)*B(:,:)
            end do
         end do

     end function la_kronecker_product_s

     ! Hermitian transpose of a matrix
     pure function la_hermitian_s(a) result(ah)
         !> Input matrix
         real(sp),intent(in) :: a(:,:)
         !> Hermitian transpose of the input matrix
         real(sp) :: ah(size(a,2,kind=ilp),size(a,1,kind=ilp))

         ah = transpose(a)

     end function la_hermitian_s

     ! Construct a square matrix whose k-th diagonal is filled with the input array
     pure function la_diag_d_from_array_k(source,k) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(dp),intent(in) :: source(:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return matrix
         real(dp) :: diag(size(source,kind=ilp) + abs(k),size(source,kind=ilp) + abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = size(source,kind=ilp)

         diag = 0.0_dp

         if (k > 0) then
            do i = 1,n
               diag(i,k + i) = source(i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i + abs(k),i) = source(i)
            end do
         else
            do i = 1,n
               diag(i,i) = source(i)
            end do
         end if

     end function la_diag_d_from_array_k

     ! Extract the main diagonal of a matrix
     pure function la_diag_d_from_matrix(A) result(diag)
         !> Input matrix. Used to define the return type and size.
         real(dp),intent(in) :: A(:,:)
         !> Return array
         real(dp) :: diag(minval(shape(A,kind=ilp)))

         !> Local variables
         integer(ilp) :: i

         do i = 1,minval(shape(A,kind=ilp))
            diag(i) = A(i,i)
         end do

     end function la_diag_d_from_matrix

     ! Extract the k-th diagonal of a matrix
     pure function la_diag_d_from_matrix_k(A,k) result(diag)
         !> Input matrix. Used to define the return type and size.
         real(dp),intent(in) :: A(:,:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return array
         real(dp) :: diag(minval(shape(A,kind=ilp)) - abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = minval(shape(A,kind=ilp)) - abs(k)

         if (k > 0) then
            do i = 1,n
               diag(i) = A(i,k + i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i) = A(i + abs(k),i)
            end do
         else
            do i = 1,n
               diag(i) = A(i,i)
            end do
         end if

     end function la_diag_d_from_matrix_k

     ! Sum of the main diagonal elements of a matrix
     pure function la_trace_d(A) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)
         !> Trace of the input matrix
         real(dp) :: res

         !> Local variables
         integer(ilp) :: i

         res = 0.0_dp
         do i = 1,minval(shape(A,kind=ilp))
            res = res + A(i,i)
         end do

     end function la_trace_d

     ! Outer product of two vectors
     pure function la_outer_product_d(u,v) result(res)
         !> First input vector
         real(dp),intent(in) :: u(:)
         !> Second input vector
         real(dp),intent(in) :: v(:)
         !> Outer product matrix
         real(dp) :: res(size(u,kind=ilp),size(v,kind=ilp))

         !> Local variables
         integer(ilp) :: col

         do col = 1,size(v,kind=ilp)
            res(:,col) = v(col)*u
         end do

     end function la_outer_product_d

     ! Cross product of two 3-dimensional vectors
     pure function la_cross_product_d(a,b) result(res)
         !> First input vector
         real(dp),intent(in) :: a(3)
         !> Second input vector
         real(dp),intent(in) :: b(3)
         !> Cross product vector
         real(dp) :: res(3)

         res(1) = a(2)*b(3) - a(3)*b(2)
         res(2) = a(3)*b(1) - a(1)*b(3)
         res(3) = a(1)*b(2) - a(2)*b(1)

     end function la_cross_product_d

     ! Kronecker product of two matrices
     pure function la_kronecker_product_d(A,B) result(C)
         !> First input matrix
         real(dp),intent(in) :: A(:,:)
         !> Second input matrix
         real(dp),intent(in) :: B(:,:)
         !> Kronecker product matrix
         real(dp) :: C(size(A,1,kind=ilp)*size(B,1,kind=ilp),size(A,2,kind=ilp)*size(B,2,kind=ilp))

         !> Local variables
         integer(ilp) :: m1,n1,maxM1,maxN1,maxM2,maxN2

         maxM1 = size(A,1,kind=ilp)
         maxN1 = size(A,2,kind=ilp)
         maxM2 = size(B,1,kind=ilp)
         maxN2 = size(B,2,kind=ilp)

         do n1 = 1,maxN1
            do m1 = 1,maxM1
               C((m1 - 1)*maxM2 + 1:m1*maxM2, (n1 - 1)*maxN2 + 1:n1*maxN2) = A(m1,n1)*B(:,:)
            end do
         end do

     end function la_kronecker_product_d

     ! Hermitian transpose of a matrix
     pure function la_hermitian_d(a) result(ah)
         !> Input matrix
         real(dp),intent(in) :: a(:,:)
         !> Hermitian transpose of the input matrix
         real(dp) :: ah(size(a,2,kind=ilp),size(a,1,kind=ilp))

         ah = transpose(a)

     end function la_hermitian_d

     ! Construct a square matrix whose k-th diagonal is filled with the input array
     pure function la_diag_q_from_array_k(source,k) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         real(qp),intent(in) :: source(:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return matrix
         real(qp) :: diag(size(source,kind=ilp) + abs(k),size(source,kind=ilp) + abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = size(source,kind=ilp)

         diag = 0.0_qp

         if (k > 0) then
            do i = 1,n
               diag(i,k + i) = source(i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i + abs(k),i) = source(i)
            end do
         else
            do i = 1,n
               diag(i,i) = source(i)
            end do
         end if

     end function la_diag_q_from_array_k

     ! Extract the main diagonal of a matrix
     pure function la_diag_q_from_matrix(A) result(diag)
         !> Input matrix. Used to define the return type and size.
         real(qp),intent(in) :: A(:,:)
         !> Return array
         real(qp) :: diag(minval(shape(A,kind=ilp)))

         !> Local variables
         integer(ilp) :: i

         do i = 1,minval(shape(A,kind=ilp))
            diag(i) = A(i,i)
         end do

     end function la_diag_q_from_matrix

     ! Extract the k-th diagonal of a matrix
     pure function la_diag_q_from_matrix_k(A,k) result(diag)
         !> Input matrix. Used to define the return type and size.
         real(qp),intent(in) :: A(:,:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return array
         real(qp) :: diag(minval(shape(A,kind=ilp)) - abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = minval(shape(A,kind=ilp)) - abs(k)

         if (k > 0) then
            do i = 1,n
               diag(i) = A(i,k + i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i) = A(i + abs(k),i)
            end do
         else
            do i = 1,n
               diag(i) = A(i,i)
            end do
         end if

     end function la_diag_q_from_matrix_k

     ! Sum of the main diagonal elements of a matrix
     pure function la_trace_q(A) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)
         !> Trace of the input matrix
         real(qp) :: res

         !> Local variables
         integer(ilp) :: i

         res = 0.0_qp
         do i = 1,minval(shape(A,kind=ilp))
            res = res + A(i,i)
         end do

     end function la_trace_q

     ! Outer product of two vectors
     pure function la_outer_product_q(u,v) result(res)
         !> First input vector
         real(qp),intent(in) :: u(:)
         !> Second input vector
         real(qp),intent(in) :: v(:)
         !> Outer product matrix
         real(qp) :: res(size(u,kind=ilp),size(v,kind=ilp))

         !> Local variables
         integer(ilp) :: col

         do col = 1,size(v,kind=ilp)
            res(:,col) = v(col)*u
         end do

     end function la_outer_product_q

     ! Cross product of two 3-dimensional vectors
     pure function la_cross_product_q(a,b) result(res)
         !> First input vector
         real(qp),intent(in) :: a(3)
         !> Second input vector
         real(qp),intent(in) :: b(3)
         !> Cross product vector
         real(qp) :: res(3)

         res(1) = a(2)*b(3) - a(3)*b(2)
         res(2) = a(3)*b(1) - a(1)*b(3)
         res(3) = a(1)*b(2) - a(2)*b(1)

     end function la_cross_product_q

     ! Kronecker product of two matrices
     pure function la_kronecker_product_q(A,B) result(C)
         !> First input matrix
         real(qp),intent(in) :: A(:,:)
         !> Second input matrix
         real(qp),intent(in) :: B(:,:)
         !> Kronecker product matrix
         real(qp) :: C(size(A,1,kind=ilp)*size(B,1,kind=ilp),size(A,2,kind=ilp)*size(B,2,kind=ilp))

         !> Local variables
         integer(ilp) :: m1,n1,maxM1,maxN1,maxM2,maxN2

         maxM1 = size(A,1,kind=ilp)
         maxN1 = size(A,2,kind=ilp)
         maxM2 = size(B,1,kind=ilp)
         maxN2 = size(B,2,kind=ilp)

         do n1 = 1,maxN1
            do m1 = 1,maxM1
               C((m1 - 1)*maxM2 + 1:m1*maxM2, (n1 - 1)*maxN2 + 1:n1*maxN2) = A(m1,n1)*B(:,:)
            end do
         end do

     end function la_kronecker_product_q

     ! Hermitian transpose of a matrix
     pure function la_hermitian_q(a) result(ah)
         !> Input matrix
         real(qp),intent(in) :: a(:,:)
         !> Hermitian transpose of the input matrix
         real(qp) :: ah(size(a,2,kind=ilp),size(a,1,kind=ilp))

         ah = transpose(a)

     end function la_hermitian_q

     ! Construct a square matrix whose k-th diagonal is filled with the input array
     pure function la_diag_c_from_array_k(source,k) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(sp),intent(in) :: source(:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return matrix
         complex(sp) :: diag(size(source,kind=ilp) + abs(k),size(source,kind=ilp) + abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = size(source,kind=ilp)

         diag = 0.0_sp

         if (k > 0) then
            do i = 1,n
               diag(i,k + i) = source(i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i + abs(k),i) = source(i)
            end do
         else
            do i = 1,n
               diag(i,i) = source(i)
            end do
         end if

     end function la_diag_c_from_array_k

     ! Extract the main diagonal of a matrix
     pure function la_diag_c_from_matrix(A) result(diag)
         !> Input matrix. Used to define the return type and size.
         complex(sp),intent(in) :: A(:,:)
         !> Return array
         complex(sp) :: diag(minval(shape(A,kind=ilp)))

         !> Local variables
         integer(ilp) :: i

         do i = 1,minval(shape(A,kind=ilp))
            diag(i) = A(i,i)
         end do

     end function la_diag_c_from_matrix

     ! Extract the k-th diagonal of a matrix
     pure function la_diag_c_from_matrix_k(A,k) result(diag)
         !> Input matrix. Used to define the return type and size.
         complex(sp),intent(in) :: A(:,:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return array
         complex(sp) :: diag(minval(shape(A,kind=ilp)) - abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = minval(shape(A,kind=ilp)) - abs(k)

         if (k > 0) then
            do i = 1,n
               diag(i) = A(i,k + i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i) = A(i + abs(k),i)
            end do
         else
            do i = 1,n
               diag(i) = A(i,i)
            end do
         end if

     end function la_diag_c_from_matrix_k

     ! Sum of the main diagonal elements of a matrix
     pure function la_trace_c(A) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)
         !> Trace of the input matrix
         complex(sp) :: res

         !> Local variables
         integer(ilp) :: i

         res = 0.0_sp
         do i = 1,minval(shape(A,kind=ilp))
            res = res + A(i,i)
         end do

     end function la_trace_c

     ! Outer product of two vectors
     pure function la_outer_product_c(u,v) result(res)
         !> First input vector
         complex(sp),intent(in) :: u(:)
         !> Second input vector
         complex(sp),intent(in) :: v(:)
         !> Outer product matrix
         complex(sp) :: res(size(u,kind=ilp),size(v,kind=ilp))

         !> Local variables
         integer(ilp) :: col

         do col = 1,size(v,kind=ilp)
            res(:,col) = v(col)*u
         end do

     end function la_outer_product_c

     ! Cross product of two 3-dimensional vectors
     pure function la_cross_product_c(a,b) result(res)
         !> First input vector
         complex(sp),intent(in) :: a(3)
         !> Second input vector
         complex(sp),intent(in) :: b(3)
         !> Cross product vector
         complex(sp) :: res(3)

         res(1) = a(2)*b(3) - a(3)*b(2)
         res(2) = a(3)*b(1) - a(1)*b(3)
         res(3) = a(1)*b(2) - a(2)*b(1)

     end function la_cross_product_c

     ! Kronecker product of two matrices
     pure function la_kronecker_product_c(A,B) result(C)
         !> First input matrix
         complex(sp),intent(in) :: A(:,:)
         !> Second input matrix
         complex(sp),intent(in) :: B(:,:)
         !> Kronecker product matrix
         complex(sp) :: C(size(A,1,kind=ilp)*size(B,1,kind=ilp),size(A,2,kind=ilp)*size(B,2,kind=ilp))

         !> Local variables
         integer(ilp) :: m1,n1,maxM1,maxN1,maxM2,maxN2

         maxM1 = size(A,1,kind=ilp)
         maxN1 = size(A,2,kind=ilp)
         maxM2 = size(B,1,kind=ilp)
         maxN2 = size(B,2,kind=ilp)

         do n1 = 1,maxN1
            do m1 = 1,maxM1
               C((m1 - 1)*maxM2 + 1:m1*maxM2, (n1 - 1)*maxN2 + 1:n1*maxN2) = A(m1,n1)*B(:,:)
            end do
         end do

     end function la_kronecker_product_c

     ! Hermitian transpose of a matrix
     pure function la_hermitian_c(a) result(ah)
         !> Input matrix
         complex(sp),intent(in) :: a(:,:)
         !> Hermitian transpose of the input matrix
         complex(sp) :: ah(size(a,2,kind=ilp),size(a,1,kind=ilp))

         ah = conjg(transpose(a))

     end function la_hermitian_c

     ! Construct a square matrix whose k-th diagonal is filled with the input array
     pure function la_diag_z_from_array_k(source,k) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(dp),intent(in) :: source(:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return matrix
         complex(dp) :: diag(size(source,kind=ilp) + abs(k),size(source,kind=ilp) + abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = size(source,kind=ilp)

         diag = 0.0_dp

         if (k > 0) then
            do i = 1,n
               diag(i,k + i) = source(i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i + abs(k),i) = source(i)
            end do
         else
            do i = 1,n
               diag(i,i) = source(i)
            end do
         end if

     end function la_diag_z_from_array_k

     ! Extract the main diagonal of a matrix
     pure function la_diag_z_from_matrix(A) result(diag)
         !> Input matrix. Used to define the return type and size.
         complex(dp),intent(in) :: A(:,:)
         !> Return array
         complex(dp) :: diag(minval(shape(A,kind=ilp)))

         !> Local variables
         integer(ilp) :: i

         do i = 1,minval(shape(A,kind=ilp))
            diag(i) = A(i,i)
         end do

     end function la_diag_z_from_matrix

     ! Extract the k-th diagonal of a matrix
     pure function la_diag_z_from_matrix_k(A,k) result(diag)
         !> Input matrix. Used to define the return type and size.
         complex(dp),intent(in) :: A(:,:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return array
         complex(dp) :: diag(minval(shape(A,kind=ilp)) - abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = minval(shape(A,kind=ilp)) - abs(k)

         if (k > 0) then
            do i = 1,n
               diag(i) = A(i,k + i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i) = A(i + abs(k),i)
            end do
         else
            do i = 1,n
               diag(i) = A(i,i)
            end do
         end if

     end function la_diag_z_from_matrix_k

     ! Sum of the main diagonal elements of a matrix
     pure function la_trace_z(A) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)
         !> Trace of the input matrix
         complex(dp) :: res

         !> Local variables
         integer(ilp) :: i

         res = 0.0_dp
         do i = 1,minval(shape(A,kind=ilp))
            res = res + A(i,i)
         end do

     end function la_trace_z

     ! Outer product of two vectors
     pure function la_outer_product_z(u,v) result(res)
         !> First input vector
         complex(dp),intent(in) :: u(:)
         !> Second input vector
         complex(dp),intent(in) :: v(:)
         !> Outer product matrix
         complex(dp) :: res(size(u,kind=ilp),size(v,kind=ilp))

         !> Local variables
         integer(ilp) :: col

         do col = 1,size(v,kind=ilp)
            res(:,col) = v(col)*u
         end do

     end function la_outer_product_z

     ! Cross product of two 3-dimensional vectors
     pure function la_cross_product_z(a,b) result(res)
         !> First input vector
         complex(dp),intent(in) :: a(3)
         !> Second input vector
         complex(dp),intent(in) :: b(3)
         !> Cross product vector
         complex(dp) :: res(3)

         res(1) = a(2)*b(3) - a(3)*b(2)
         res(2) = a(3)*b(1) - a(1)*b(3)
         res(3) = a(1)*b(2) - a(2)*b(1)

     end function la_cross_product_z

     ! Kronecker product of two matrices
     pure function la_kronecker_product_z(A,B) result(C)
         !> First input matrix
         complex(dp),intent(in) :: A(:,:)
         !> Second input matrix
         complex(dp),intent(in) :: B(:,:)
         !> Kronecker product matrix
         complex(dp) :: C(size(A,1,kind=ilp)*size(B,1,kind=ilp),size(A,2,kind=ilp)*size(B,2,kind=ilp))

         !> Local variables
         integer(ilp) :: m1,n1,maxM1,maxN1,maxM2,maxN2

         maxM1 = size(A,1,kind=ilp)
         maxN1 = size(A,2,kind=ilp)
         maxM2 = size(B,1,kind=ilp)
         maxN2 = size(B,2,kind=ilp)

         do n1 = 1,maxN1
            do m1 = 1,maxM1
               C((m1 - 1)*maxM2 + 1:m1*maxM2, (n1 - 1)*maxN2 + 1:n1*maxN2) = A(m1,n1)*B(:,:)
            end do
         end do

     end function la_kronecker_product_z

     ! Hermitian transpose of a matrix
     pure function la_hermitian_z(a) result(ah)
         !> Input matrix
         complex(dp),intent(in) :: a(:,:)
         !> Hermitian transpose of the input matrix
         complex(dp) :: ah(size(a,2,kind=ilp),size(a,1,kind=ilp))

         ah = conjg(transpose(a))

     end function la_hermitian_z

     ! Construct a square matrix whose k-th diagonal is filled with the input array
     pure function la_diag_w_from_array_k(source,k) result(diag)
         !> Array of diagonal values. Used to define the return type and the matrix size.
         complex(qp),intent(in) :: source(:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return matrix
         complex(qp) :: diag(size(source,kind=ilp) + abs(k),size(source,kind=ilp) + abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = size(source,kind=ilp)

         diag = 0.0_qp

         if (k > 0) then
            do i = 1,n
               diag(i,k + i) = source(i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i + abs(k),i) = source(i)
            end do
         else
            do i = 1,n
               diag(i,i) = source(i)
            end do
         end if

     end function la_diag_w_from_array_k

     ! Extract the main diagonal of a matrix
     pure function la_diag_w_from_matrix(A) result(diag)
         !> Input matrix. Used to define the return type and size.
         complex(qp),intent(in) :: A(:,:)
         !> Return array
         complex(qp) :: diag(minval(shape(A,kind=ilp)))

         !> Local variables
         integer(ilp) :: i

         do i = 1,minval(shape(A,kind=ilp))
            diag(i) = A(i,i)
         end do

     end function la_diag_w_from_matrix

     ! Extract the k-th diagonal of a matrix
     pure function la_diag_w_from_matrix_k(A,k) result(diag)
         !> Input matrix. Used to define the return type and size.
         complex(qp),intent(in) :: A(:,:)
         !> Index of the diagonal: 0 is the main diagonal, k>0 the k-th superdiagonal, k<0 the k-th subdiagonal.
         integer(ilp),intent(in) :: k
         !> Return array
         complex(qp) :: diag(minval(shape(A,kind=ilp)) - abs(k))

         !> Local variables
         integer(ilp) :: i,n

         n = minval(shape(A,kind=ilp)) - abs(k)

         if (k > 0) then
            do i = 1,n
               diag(i) = A(i,k + i)
            end do
         elseif (k < 0) then
            do i = 1,n
               diag(i) = A(i + abs(k),i)
            end do
         else
            do i = 1,n
               diag(i) = A(i,i)
            end do
         end if

     end function la_diag_w_from_matrix_k

     ! Sum of the main diagonal elements of a matrix
     pure function la_trace_w(A) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)
         !> Trace of the input matrix
         complex(qp) :: res

         !> Local variables
         integer(ilp) :: i

         res = 0.0_qp
         do i = 1,minval(shape(A,kind=ilp))
            res = res + A(i,i)
         end do

     end function la_trace_w

     ! Outer product of two vectors
     pure function la_outer_product_w(u,v) result(res)
         !> First input vector
         complex(qp),intent(in) :: u(:)
         !> Second input vector
         complex(qp),intent(in) :: v(:)
         !> Outer product matrix
         complex(qp) :: res(size(u,kind=ilp),size(v,kind=ilp))

         !> Local variables
         integer(ilp) :: col

         do col = 1,size(v,kind=ilp)
            res(:,col) = v(col)*u
         end do

     end function la_outer_product_w

     ! Cross product of two 3-dimensional vectors
     pure function la_cross_product_w(a,b) result(res)
         !> First input vector
         complex(qp),intent(in) :: a(3)
         !> Second input vector
         complex(qp),intent(in) :: b(3)
         !> Cross product vector
         complex(qp) :: res(3)

         res(1) = a(2)*b(3) - a(3)*b(2)
         res(2) = a(3)*b(1) - a(1)*b(3)
         res(3) = a(1)*b(2) - a(2)*b(1)

     end function la_cross_product_w

     ! Kronecker product of two matrices
     pure function la_kronecker_product_w(A,B) result(C)
         !> First input matrix
         complex(qp),intent(in) :: A(:,:)
         !> Second input matrix
         complex(qp),intent(in) :: B(:,:)
         !> Kronecker product matrix
         complex(qp) :: C(size(A,1,kind=ilp)*size(B,1,kind=ilp),size(A,2,kind=ilp)*size(B,2,kind=ilp))

         !> Local variables
         integer(ilp) :: m1,n1,maxM1,maxN1,maxM2,maxN2

         maxM1 = size(A,1,kind=ilp)
         maxN1 = size(A,2,kind=ilp)
         maxM2 = size(B,1,kind=ilp)
         maxN2 = size(B,2,kind=ilp)

         do n1 = 1,maxN1
            do m1 = 1,maxM1
               C((m1 - 1)*maxM2 + 1:m1*maxM2, (n1 - 1)*maxN2 + 1:n1*maxN2) = A(m1,n1)*B(:,:)
            end do
         end do

     end function la_kronecker_product_w

     ! Hermitian transpose of a matrix
     pure function la_hermitian_w(a) result(ah)
         !> Input matrix
         complex(qp),intent(in) :: a(:,:)
         !> Hermitian transpose of the input matrix
         complex(qp) :: ah(size(a,2,kind=ilp),size(a,1,kind=ilp))

         ah = conjg(transpose(a))

     end function la_hermitian_w

end module la_eye
