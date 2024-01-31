module stdlib_linalg_state
     use stdlib_linalg_constants
     implicit none(type, external)
     public

     !> State return types
     integer(ilp), parameter :: LINALG_SUCCESS = 0_ilp

     !> Use fixed-size character storage for performance
     integer(ilp), parameter, private :: MSG_LENGTH  = 256_ilp
     integer(ilp), parameter, private :: NAME_LENGTH = 32_ilp

     !> `stdlib_state` defines a state return type for a
     !> linear algebra routine
     type, public :: stdlib_state

         !> The current exit state
         integer(ilp) :: state = LINALG_SUCCESS

         !> Message associated to the current state
         character(len=MSG_LENGTH) :: message = repeat(' ',MSG_LENGTH)

         !> Location of the state change
         character(len=NAME_LENGTH) :: where_at = repeat(' ',NAME_LENGTH)


         contains

            !> State properties
            procedure :: ok    => state_is_ok
            procedure :: error => state_is_error

     end type stdlib_state


     contains

     !> Check if the current state is successful
     elemental logical(lk) function state_is_ok(this)
        class(stdlib_state), intent(in) :: this
        state_is_ok = this%state==LINALG_SUCCESS
     end function state_is_ok

     !> Check if the current state is an error state
     elemental logical(lk) function state_is_error(this)
        class(stdlib_state), intent(in) :: this
        state_is_error = this%state/=LINALG_SUCCESS
     end function state_is_error

end module stdlib_linalg_state
