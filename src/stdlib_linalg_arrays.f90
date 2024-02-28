!> Get stride/contiguous attribute from an array using ISO_Fortran_binding.h
module stdlib_linalg_arrays
     use stdlib_linalg_constants,only:ilp,lk
     use iso_c_binding
     use iso_fortran_env, only: output_unit
     implicit none(type,external)
     public

     integer, parameter :: CFI_RANK = c_int8_t
     integer, parameter :: CFI_FLAG = c_int
     integer, parameter :: CFI_TYPE = c_int16_t
     integer, parameter :: CFI_SIZE = c_intptr_t

     ! CFI_attribute_t is a typedef name for a standard integer type capable
     ! of representing the values of the attribute codes.
     integer, parameter :: CFI_ATTR = c_int8_t

     integer(CFI_RANK), parameter, public :: MAX_FORTRAN_RANK = 15_CFI_RANK
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INTEGER = 1
     integer(CFI_RANK), parameter, public :: CFI_TYPES_LOGICAL = 2
     integer(CFI_RANK), parameter, public :: CFI_TYPES_REAL = 3
     integer(CFI_RANK), parameter, public :: CFI_TYPES_COMPLEX = 4
     integer(CFI_RANK), parameter, public :: CFI_TYPES_CHARACTER = 5
     integer(CFI_RANK), parameter, public :: CFI_TYPES_STRUCT = 6
     integer(CFI_RANK), parameter, public :: CFI_TYPES_CPTR = 7
     integer(CFI_RANK), parameter, public :: CFI_TYPES_CFUNPTR = 8
     integer(CFI_RANK), parameter, public :: CFI_TYPES_OTHER = 9
     integer(CFI_RANK), parameter, public :: CFI_TYPES_SIGNED_CHAR = 10
     integer(CFI_RANK), parameter, public :: CFI_TYPES_SHORT = 11
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT = 12
     integer(CFI_RANK), parameter, public :: CFI_TYPES_LONG = 13
     integer(CFI_RANK), parameter, public :: CFI_TYPES_LONG_LONG = 14
     integer(CFI_RANK), parameter, public :: CFI_TYPES_SIZE_T = 15
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT8_T = 16
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT16_T = 17
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT32_T = 18
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT64_T = 19
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_LEAST8_T = 20
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_LEAST16_T = 21
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_LEAST32_T = 22
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_LEAST64_T = 23
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_FAST8_T = 24
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_FAST16_T = 25
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_FAST32_T = 26
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INT_FAST64_T = 27
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INTMAX_T = 28
     integer(CFI_RANK), parameter, public :: CFI_TYPES_INTPTR_T = 29
     integer(CFI_RANK), parameter, public :: CFI_TYPES_PTRDIFF_T = 30
     integer(CFI_RANK), parameter, public :: CFI_TYPES_BOOL = 31
     integer(CFI_RANK), parameter, public :: CFI_TYPES_FLOAT = 32
     integer(CFI_RANK), parameter, public :: CFI_TYPES_DOUBLE = 33
     integer(CFI_RANK), parameter, public :: CFI_TYPES_FLOAT_COMPLEX = 34
     integer(CFI_RANK), parameter, public :: CFI_TYPES_DOUBLE_COMPLEX = 35

     !> Reproduces CFI_dim_t
     type, public, bind(C) :: array_dimension
        integer(CFI_SIZE) :: lower_bound  = 0_CFI_SIZE
        integer(CFI_SIZE) :: extent       = 0_CFI_SIZE
        integer(CFI_SIZE) :: stride_bytes = 0_CFI_SIZE
     end type array_dimension

     !> A Derived type containing array information
     type, public, bind(C) :: array_descriptor
        type(c_ptr)        :: base_addr = c_null_ptr
        integer(CFI_SIZE)  :: elem_bytes = 0_CFI_SIZE
        integer(CFI_FLAG)  :: version = 0_CFI_FLAG
        integer(CFI_RANK)  :: rank = 0_CFI_RANK
        integer(CFI_TYPE)  :: type = CFI_TYPES_OTHER
        integer(CFI_ATTR)  :: attribute = 0_CFI_ATTR
        type(array_dimension) :: dim(MAX_FORTRAN_RANK)

     end type array_descriptor

     interface array_descriptor
        module procedure get_array_descriptor
     end interface array_descriptor

     !> C Interface
     interface
         type(array_descriptor) function CFI_to_Fortran(descr) bind(C,name="CFI_to_Fortran")
            import array_descriptor
            type(*), dimension(..), intent(inout) :: descr
         end function CFI_to_fortran
     end interface


     contains

     !> Return array stride information
     integer(CFI_SIZE) function stride(array)
         type(*), dimension(:), intent(inout), target :: array
         type(array_descriptor) :: CFI
         real :: rstride
         CFI = array_descriptor(array)

         !> Number of elements between two elements of dimension i.
         !> This gives the exact number of allocated elements in all previous dimensions
         rstride  = real(CFI%dim(1)%stride_bytes) / real(CFI%elem_bytes)
         stride   = nint(rstride,kind=CFI_SIZE)

         !> Non integer stride: there is some padding
         if (.not.rstride-stride<0.00001 .or. CFI%rank/=1) stride=-1_CFI_SIZE

     end function stride

     !> Return leading dimension of a (possibly-sliced) 2d matrix
     integer(CFI_SIZE) function leading_dim(matrix)
         type(*), dimension(:,:), intent(inout), target :: matrix
         type(array_descriptor) :: CFI
         real :: rstride
         CFI = array_descriptor(matrix)

         !> Number of elements between two elements of dimension i.
         !> This gives the exact number of allocated elements in all previous dimensions
         rstride       = real(CFI%dim(2)%stride_bytes) / real(CFI%elem_bytes)
         leading_dim   = nint(rstride,kind=CFI_SIZE)

         !> Non integer stride: there is some padding
         if (.not.rstride-leading_dim<0.00001 .or. CFI%rank/=2 .or. leading_dim<0) leading_dim=0_CFI_SIZE

     end function leading_dim

     !> Return CFI descriptor corresponding to a Fortran variable
     type(array_descriptor) function get_array_descriptor(variable)
        type(*), dimension(..), intent(inout), target :: variable
        get_array_descriptor = CFI_to_Fortran(variable)
     end function get_array_descriptor

     !> print information
     subroutine CFI_print(this,unit)
        type(array_descriptor), intent(inout), target :: this
        integer, optional, intent(in) :: unit

        integer :: useUnit,i
        integer(CFI_SIZE), allocatable :: dim_strides(:)

        if (present(unit)) then
            useUnit = unit
        else
            useUnit = output_unit
        end if

        write(useUnit,1) c_loc(this%base_addr)
        write(useUnit,2) this%elem_bytes
        write(useUnit,3) this%version
        write(useUnit,4) this%type
        write(useUnit,5) this%rank

        !if (this%rank>0) dim_strides = CFI_strides(this)

        do i=1,this%rank

!            select case (dim_strides(i))
!               case (1)
!                   write(useUnit,7) i,this%dim(i)%lower_bound+1, & ! C->Fortran indexing
!                                      this%dim(i)%lower_bound+this%dim(i)%extent
!               case (-1)
                   ! Not an integer number of elements
                   write(useUnit,7) i,this%dim(i)%lower_bound+1, & ! C->Fortran indexing
                                      this%dim(i)%lower_bound+this%dim(i)%extent, &
                                      nint(this%dim(i)%stride_bytes/real(this%elem_bytes))
!               case default
!                   ! The stride is exactly an integer number of elements
!                   write(useUnit,6) i,this%dim(i)%lower_bound+1, & ! C->Fortran indexing
!                                      this%dim(i)%lower_bound+dim_strides(i)*this%dim(i)%extent, &
!                                      dim_strides(i)
!
!            end select

        end do

        1 format('base address  : ',z0)
        2 format('total length  : ',i0)
        3 format('CFI version   : ',i0)
        4 format('Variable type : ',i0,:,1x,'(',a,')')
        5 format('Rank          : ',i0)
        6 format('  ',i2,') [',i0,':',i0,':',i0,']') ! Use fortran syntax: stride is last
        7 format('  ',i2,') [',i0,':',i0,']',:,', padded every ',i0,' elements')

     end subroutine CFI_print


end module stdlib_linalg_arrays

