!> Matrix property checks
module la_matrix_property_checks
     use la_constants
     use la_state_type
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> @brief Check whether a matrix is square.
     !!
     !! This function returns `.true.` if the input matrix has as many rows as columns.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !!
     !! @return `.true.` if \f$ m = n \f$, `.false.` otherwise.
     !!
     public :: is_square

     !> @brief Check whether a matrix is diagonal.
     !!
     !! This function returns `.true.` if every off-diagonal entry of the input matrix is zero.
     !! The matrix need not be square.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !!
     !! @return `.true.` if \f$ A_{ij} = 0 \f$ for all \f$ i \neq j \f$, `.false.` otherwise.
     !!
     public :: is_diagonal

     !> @brief Check whether a matrix is symmetric.
     !!
     !! This function returns `.true.` if the input matrix equals its own transpose.
     !! A non-square matrix is never symmetric.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !!
     !! @return `.true.` if \f$ A = A^T \f$, `.false.` otherwise.
     !!
     public :: is_symmetric

     !> @brief Check whether a matrix is skew-symmetric.
     !!
     !! This function returns `.true.` if the input matrix equals the negative of its own transpose.
     !! A non-square matrix is never skew-symmetric.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !!
     !! @return `.true.` if \f$ A = -A^T \f$, `.false.` otherwise.
     !!
     public :: is_skew_symmetric

     !> @brief Check whether a matrix is Hermitian.
     !!
     !! This function returns `.true.` if the input matrix equals its own conjugate transpose.
     !! For a real matrix this is the same test as [`is_symmetric`](@ref la_matrix_property_checks::is_symmetric).
     !! A non-square matrix is never Hermitian.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !!
     !! @return `.true.` if \f$ A = A^H \f$, `.false.` otherwise.
     !!
     public :: is_hermitian

     !> @brief Check whether a matrix is triangular.
     !!
     !! This function returns `.true.` if every entry of the input matrix below (`uplo = 'U'`) or
     !! above (`uplo = 'L'`) the main diagonal is zero. The matrix need not be square.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !! @param[in] uplo The triangle to test: `'U'` for upper, `'L'` for lower.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return `.true.` if `A` is triangular of the requested type, `.false.` otherwise.
     !!
     public :: is_triangular

     !> @brief Check whether a matrix is Hessenberg.
     !!
     !! This function returns `.true.` if every entry of the input matrix more than one row below
     !! (`uplo = 'U'`) or more than one row above (`uplo = 'L'`) the main diagonal is zero.
     !! The matrix need not be square.
     !!
     !! @param[in] A The input matrix of size \f$ [m,n] \f$.
     !! @param[in] uplo The Hessenberg form to test: `'U'` for upper, `'L'` for lower.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return `.true.` if `A` is Hessenberg of the requested type, `.false.` otherwise.
     !!
     public :: is_hessenberg

     ! Exact inequality of two values. It is written as a relational test on the magnitude of the
     ! difference, which the compiler does not flag as a floating-point equality comparison, and
     ! which reports NaN operands as different exactly as `/=` does.
     interface differ
        module procedure la_differ_s
        module procedure la_differ_d
        module procedure la_differ_q
        module procedure la_differ_c
        module procedure la_differ_z
        module procedure la_differ_w
     end interface differ

     interface is_square
        module procedure la_is_square_s
        module procedure la_is_square_d
        module procedure la_is_square_q
        module procedure la_is_square_c
        module procedure la_is_square_z
        module procedure la_is_square_w
     end interface is_square

     interface is_diagonal
        module procedure la_is_diagonal_s
        module procedure la_is_diagonal_d
        module procedure la_is_diagonal_q
        module procedure la_is_diagonal_c
        module procedure la_is_diagonal_z
        module procedure la_is_diagonal_w
     end interface is_diagonal

     interface is_symmetric
        module procedure la_is_symmetric_s
        module procedure la_is_symmetric_d
        module procedure la_is_symmetric_q
        module procedure la_is_symmetric_c
        module procedure la_is_symmetric_z
        module procedure la_is_symmetric_w
     end interface is_symmetric

     interface is_skew_symmetric
        module procedure la_is_skew_symmetric_s
        module procedure la_is_skew_symmetric_d
        module procedure la_is_skew_symmetric_q
        module procedure la_is_skew_symmetric_c
        module procedure la_is_skew_symmetric_z
        module procedure la_is_skew_symmetric_w
     end interface is_skew_symmetric

     interface is_hermitian
        module procedure la_is_hermitian_s
        module procedure la_is_hermitian_d
        module procedure la_is_hermitian_q
        module procedure la_is_hermitian_c
        module procedure la_is_hermitian_z
        module procedure la_is_hermitian_w
     end interface is_hermitian

     interface is_triangular
        module procedure la_is_triangular_s
        module procedure la_is_triangular_d
        module procedure la_is_triangular_q
        module procedure la_is_triangular_c
        module procedure la_is_triangular_z
        module procedure la_is_triangular_w
        module procedure la_is_triangular_s_errhandle
        module procedure la_is_triangular_d_errhandle
        module procedure la_is_triangular_q_errhandle
        module procedure la_is_triangular_c_errhandle
        module procedure la_is_triangular_z_errhandle
        module procedure la_is_triangular_w_errhandle
     end interface is_triangular

     interface is_hessenberg
        module procedure la_is_hessenberg_s
        module procedure la_is_hessenberg_d
        module procedure la_is_hessenberg_q
        module procedure la_is_hessenberg_c
        module procedure la_is_hessenberg_z
        module procedure la_is_hessenberg_w
        module procedure la_is_hessenberg_s_errhandle
        module procedure la_is_hessenberg_d_errhandle
        module procedure la_is_hessenberg_q_errhandle
        module procedure la_is_hessenberg_c_errhandle
        module procedure la_is_hessenberg_z_errhandle
        module procedure la_is_hessenberg_w_errhandle
     end interface is_hessenberg

     contains

     pure elemental logical(lk) function la_differ_s(x,y) result(differs)
         real(sp),intent(in) :: x,y

         differs = .not. abs(x - y) <= 0.0_sp

     end function la_differ_s

     !> Check whether a matrix is square.
     pure logical(lk) function la_is_square_s(A) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)

         res = size(A,1,kind=ilp) == size(A,2,kind=ilp)

     end function la_is_square_s

     !> Check whether a matrix is diagonal.
     pure logical(lk) function la_is_diagonal_s(A) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)

         real(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)

         do j = 1,n
            o = min(j - 1,m)
            do i = 1,o
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
            do i = o + 2,m
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_diagonal_s

     !> Check whether a matrix is symmetric.
     pure logical(lk) function la_is_symmetric_s(A) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j - 1
               if (differ(A(i,j),A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_symmetric_s

     !> Check whether a matrix is skew-symmetric.
     pure logical(lk) function la_is_skew_symmetric_s(A) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),-A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_skew_symmetric_s

     pure elemental logical(lk) function la_differ_d(x,y) result(differs)
         real(dp),intent(in) :: x,y

         differs = .not. abs(x - y) <= 0.0_dp

     end function la_differ_d

     !> Check whether a matrix is square.
     pure logical(lk) function la_is_square_d(A) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)

         res = size(A,1,kind=ilp) == size(A,2,kind=ilp)

     end function la_is_square_d

     !> Check whether a matrix is diagonal.
     pure logical(lk) function la_is_diagonal_d(A) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)

         real(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)

         do j = 1,n
            o = min(j - 1,m)
            do i = 1,o
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
            do i = o + 2,m
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_diagonal_d

     !> Check whether a matrix is symmetric.
     pure logical(lk) function la_is_symmetric_d(A) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j - 1
               if (differ(A(i,j),A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_symmetric_d

     !> Check whether a matrix is skew-symmetric.
     pure logical(lk) function la_is_skew_symmetric_d(A) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),-A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_skew_symmetric_d

     pure elemental logical(lk) function la_differ_q(x,y) result(differs)
         real(qp),intent(in) :: x,y

         differs = .not. abs(x - y) <= 0.0_qp

     end function la_differ_q

     !> Check whether a matrix is square.
     pure logical(lk) function la_is_square_q(A) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)

         res = size(A,1,kind=ilp) == size(A,2,kind=ilp)

     end function la_is_square_q

     !> Check whether a matrix is diagonal.
     pure logical(lk) function la_is_diagonal_q(A) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)

         real(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)

         do j = 1,n
            o = min(j - 1,m)
            do i = 1,o
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
            do i = o + 2,m
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_diagonal_q

     !> Check whether a matrix is symmetric.
     pure logical(lk) function la_is_symmetric_q(A) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j - 1
               if (differ(A(i,j),A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_symmetric_q

     !> Check whether a matrix is skew-symmetric.
     pure logical(lk) function la_is_skew_symmetric_q(A) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),-A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_skew_symmetric_q

     pure elemental logical(lk) function la_differ_c(x,y) result(differs)
         complex(sp),intent(in) :: x,y

         differs = .not. abs(x - y) <= 0.0_sp

     end function la_differ_c

     !> Check whether a matrix is square.
     pure logical(lk) function la_is_square_c(A) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)

         res = size(A,1,kind=ilp) == size(A,2,kind=ilp)

     end function la_is_square_c

     !> Check whether a matrix is diagonal.
     pure logical(lk) function la_is_diagonal_c(A) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)

         complex(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)

         do j = 1,n
            o = min(j - 1,m)
            do i = 1,o
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
            do i = o + 2,m
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_diagonal_c

     !> Check whether a matrix is symmetric.
     pure logical(lk) function la_is_symmetric_c(A) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j - 1
               if (differ(A(i,j),A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_symmetric_c

     !> Check whether a matrix is skew-symmetric.
     pure logical(lk) function la_is_skew_symmetric_c(A) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),-A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_skew_symmetric_c

     pure elemental logical(lk) function la_differ_z(x,y) result(differs)
         complex(dp),intent(in) :: x,y

         differs = .not. abs(x - y) <= 0.0_dp

     end function la_differ_z

     !> Check whether a matrix is square.
     pure logical(lk) function la_is_square_z(A) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)

         res = size(A,1,kind=ilp) == size(A,2,kind=ilp)

     end function la_is_square_z

     !> Check whether a matrix is diagonal.
     pure logical(lk) function la_is_diagonal_z(A) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)

         complex(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)

         do j = 1,n
            o = min(j - 1,m)
            do i = 1,o
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
            do i = o + 2,m
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_diagonal_z

     !> Check whether a matrix is symmetric.
     pure logical(lk) function la_is_symmetric_z(A) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j - 1
               if (differ(A(i,j),A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_symmetric_z

     !> Check whether a matrix is skew-symmetric.
     pure logical(lk) function la_is_skew_symmetric_z(A) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),-A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_skew_symmetric_z

     pure elemental logical(lk) function la_differ_w(x,y) result(differs)
         complex(qp),intent(in) :: x,y

         differs = .not. abs(x - y) <= 0.0_qp

     end function la_differ_w

     !> Check whether a matrix is square.
     pure logical(lk) function la_is_square_w(A) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)

         res = size(A,1,kind=ilp) == size(A,2,kind=ilp)

     end function la_is_square_w

     !> Check whether a matrix is diagonal.
     pure logical(lk) function la_is_diagonal_w(A) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)

         complex(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)

         do j = 1,n
            o = min(j - 1,m)
            do i = 1,o
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
            do i = o + 2,m
               if (differ(A(i,j),zero)) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_diagonal_w

     !> Check whether a matrix is symmetric.
     pure logical(lk) function la_is_symmetric_w(A) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j - 1
               if (differ(A(i,j),A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_symmetric_w

     !> Check whether a matrix is skew-symmetric.
     pure logical(lk) function la_is_skew_symmetric_w(A) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),-A(j,i))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_skew_symmetric_w

     !> Check whether a real matrix is Hermitian: symmetry and Hermiticity coincide.
     pure logical(lk) function la_is_hermitian_s(A) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)

         res = is_symmetric(A)

     end function la_is_hermitian_s

     !> Check whether a real matrix is Hermitian: symmetry and Hermiticity coincide.
     pure logical(lk) function la_is_hermitian_d(A) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)

         res = is_symmetric(A)

     end function la_is_hermitian_d

     !> Check whether a real matrix is Hermitian: symmetry and Hermiticity coincide.
     pure logical(lk) function la_is_hermitian_q(A) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)

         res = is_symmetric(A)

     end function la_is_hermitian_q

     !> Check whether a complex matrix is Hermitian.
     pure logical(lk) function la_is_hermitian_c(A) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),conjg(A(j,i)))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_hermitian_c

     !> Check whether a complex matrix is Hermitian.
     pure logical(lk) function la_is_hermitian_z(A) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),conjg(A(j,i)))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_hermitian_z

     !> Check whether a complex matrix is Hermitian.
     pure logical(lk) function la_is_hermitian_w(A) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)

         integer(ilp) :: n,i,j

         if (.not. is_square(A)) then
            res = .false._lk
            return
         end if

         n = size(A,1,kind=ilp)

         do j = 1,n
            do i = 1,j
               if (differ(A(i,j),conjg(A(j,i)))) then
                  res = .false._lk
                  return
               end if
            end do
         end do

         res = .true._lk

     end function la_is_hermitian_w

     !> Check whether a matrix is upper or lower triangular.
     pure logical(lk) function la_is_triangular_s(A,uplo) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         real(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_triangular_s

     !> Check whether a matrix is in upper or lower Hessenberg form.
     pure logical(lk) function la_is_hessenberg_s(A,uplo) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         real(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_hessenberg_s

     !> Check whether a matrix is upper or lower triangular.
     pure logical(lk) function la_is_triangular_d(A,uplo) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         real(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_triangular_d

     !> Check whether a matrix is in upper or lower Hessenberg form.
     pure logical(lk) function la_is_hessenberg_d(A,uplo) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         real(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_hessenberg_d

     !> Check whether a matrix is upper or lower triangular.
     pure logical(lk) function la_is_triangular_q(A,uplo) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         real(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_triangular_q

     !> Check whether a matrix is in upper or lower Hessenberg form.
     pure logical(lk) function la_is_hessenberg_q(A,uplo) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         real(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_hessenberg_q

     !> Check whether a matrix is upper or lower triangular.
     pure logical(lk) function la_is_triangular_c(A,uplo) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         complex(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_triangular_c

     !> Check whether a matrix is in upper or lower Hessenberg form.
     pure logical(lk) function la_is_hessenberg_c(A,uplo) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         complex(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_hessenberg_c

     !> Check whether a matrix is upper or lower triangular.
     pure logical(lk) function la_is_triangular_z(A,uplo) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         complex(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_triangular_z

     !> Check whether a matrix is in upper or lower Hessenberg form.
     pure logical(lk) function la_is_hessenberg_z(A,uplo) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         complex(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_hessenberg_z

     !> Check whether a matrix is upper or lower triangular.
     pure logical(lk) function la_is_triangular_w(A,uplo) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         complex(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_triangular_w

     !> Check whether a matrix is in upper or lower Hessenberg form.
     pure logical(lk) function la_is_hessenberg_w(A,uplo) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo

         complex(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle()

     end function la_is_hessenberg_w

     !> Check whether a matrix is upper or lower triangular.
     logical(lk) function la_is_triangular_s_errhandle(A,uplo,err) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         real(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_triangular_s_errhandle

     !> Check whether a matrix is in upper or lower Hessenberg form.
     logical(lk) function la_is_hessenberg_s_errhandle(A,uplo,err) result(res)
         !> Input matrix
         real(sp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         real(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_hessenberg_s_errhandle

     !> Check whether a matrix is upper or lower triangular.
     logical(lk) function la_is_triangular_d_errhandle(A,uplo,err) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         real(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_triangular_d_errhandle

     !> Check whether a matrix is in upper or lower Hessenberg form.
     logical(lk) function la_is_hessenberg_d_errhandle(A,uplo,err) result(res)
         !> Input matrix
         real(dp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         real(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_hessenberg_d_errhandle

     !> Check whether a matrix is upper or lower triangular.
     logical(lk) function la_is_triangular_q_errhandle(A,uplo,err) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         real(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_triangular_q_errhandle

     !> Check whether a matrix is in upper or lower Hessenberg form.
     logical(lk) function la_is_hessenberg_q_errhandle(A,uplo,err) result(res)
         !> Input matrix
         real(qp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         real(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_hessenberg_q_errhandle

     !> Check whether a matrix is upper or lower triangular.
     logical(lk) function la_is_triangular_c_errhandle(A,uplo,err) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         complex(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_triangular_c_errhandle

     !> Check whether a matrix is in upper or lower Hessenberg form.
     logical(lk) function la_is_hessenberg_c_errhandle(A,uplo,err) result(res)
         !> Input matrix
         complex(sp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         complex(sp),parameter :: zero = 0.0_sp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_hessenberg_c_errhandle

     !> Check whether a matrix is upper or lower triangular.
     logical(lk) function la_is_triangular_z_errhandle(A,uplo,err) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         complex(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_triangular_z_errhandle

     !> Check whether a matrix is in upper or lower Hessenberg form.
     logical(lk) function la_is_hessenberg_z_errhandle(A,uplo,err) result(res)
         !> Input matrix
         complex(dp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         complex(dp),parameter :: zero = 0.0_dp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_hessenberg_z_errhandle

     !> Check whether a matrix is upper or lower triangular.
     logical(lk) function la_is_triangular_w_errhandle(A,uplo,err) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)
         !> Triangle to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         complex(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_triangular'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = o + 2,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 1,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid triangle selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_triangular_w_errhandle

     !> Check whether a matrix is in upper or lower Hessenberg form.
     logical(lk) function la_is_hessenberg_w_errhandle(A,uplo,err) result(res)
         !> Input matrix
         complex(qp),intent(in) :: A(:,:)
         !> Hessenberg form to be tested: 'U' for upper, 'L' for lower
         character,intent(in) :: uplo
         !> [optional] state return flag. On error if not requested, the code will stop
         type(la_state),intent(out) :: err

         complex(qp),parameter :: zero = 0.0_qp
         integer(ilp) :: m,n,o,i,j
         type(la_state) :: err0
         character(*),parameter :: this = 'is_hessenberg'

         m = size(A,1,kind=ilp)
         n = size(A,2,kind=ilp)
         res = .false._lk

         select case (uplo)
            case ('u','U')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = o + 4,m
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case ('l','L')
               do j = 1,n
                  o = min(j - 2,m)
                  do i = 1,o
                     if (differ(A(i,j),zero)) goto 1
                  end do
               end do
            case default
               err0 = la_state(this,LINALG_VALUE_ERROR,'invalid Hessenberg selector: uplo=',uplo)
               goto 1
         end select

         res = .true._lk

1        call err0%handle(err)

     end function la_is_hessenberg_w_errhandle

end module la_matrix_property_checks
