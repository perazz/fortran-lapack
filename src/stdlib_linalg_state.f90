module stdlib_linalg_state
     use stdlib_linalg_constants
     implicit none(type, external)
     public

     !> Public interfaces
     public :: operator(==),operator(/=)
     public :: operator(<),operator(<=)
     public :: operator(>),operator(>=)

     !> State return types
     integer(ilp), parameter :: LINALG_SUCCESS = 0_ilp

     !> Use fixed-size character storage for performance
     integer(ilp), parameter, private :: MSG_LENGTH  = 256_ilp
     integer(ilp), parameter, private :: NAME_LENGTH = 32_ilp

     !> `stdlib_state` defines a state return type for a
     !> linear algebra routine. State contains a status flag, a comment, and a
     !> procedure specifier that can be used to mark where the error happened
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

     !> Comparison operators
     interface operator(==)
         module procedure state_eq_flag
         module procedure flag_eq_state
     end interface
     interface operator(/=)
         module procedure state_neq_flag
         module procedure flag_neq_state
     end interface
     interface operator(<)
         module procedure state_lt_flag
         module procedure flag_lt_state
     end interface
     interface operator(<=)
         module procedure state_le_flag
         module procedure flag_le_state
     end interface
     interface operator(>)
         module procedure state_gt_flag
         module procedure flag_gt_state
     end interface
     interface operator(>=)
         module procedure state_ge_flag
         module procedure flag_ge_state
     end interface


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

     !> Compare an error flag with an integer
     elemental logical function state_eq_flag(err,flag)
        type(stdlib_state), intent(in) :: err
        integer, intent(in) :: flag
        state_eq_flag = err%ierr == flag
     end function state_eq_flag
     elemental logical function flag_eq_state(flag,err)
        integer, intent(in) :: flag
        type(stdlib_state), intent(in) :: err
        flag_eq_state = err%ierr == flag
     end function flag_eq_state
     elemental logical function state_neq_flag(err,flag)
        type(stdlib_state), intent(in) :: err
        integer, intent(in) :: flag
        state_neq_flag = .not.state_eq_flag(err,flag)
     end function state_neq_flag
     elemental logical function flag_neq_state(flag,err)
        integer, intent(in) :: flag
        type(stdlib_state), intent(in) :: err
        state_neq_flag = .not.state_eq_flag(err,flag)
     end function flag_neq_state
     elemental logical function state_lt_flag(err,flag)
        type(stdlib_state), intent(in) :: err
        integer, intent(in) :: flag
        state_lt_flag = err%ierr < flag
     end function state_lt_flag
     elemental logical function state_le_flag(err,flag)
        type(stdlib_state), intent(in) :: err
        integer, intent(in) :: flag
        state_le_flag = err%ierr <= flag
     end function state_le_flag
     elemental logical function flag_lt_state(flag,err)
        integer, intent(in) :: flag
        type(stdlib_state), intent(in) :: err
        flag_lt_state = err%ierr < flag
     end function flag_lt_state
     elemental logical function flag_le_state(flag,err)
        integer, intent(in) :: flag
        type(stdlib_state), intent(in) :: err
        flag_le_state = err%ierr <= flag
     end function flag_le_state
     elemental logical function state_gt_flag(err,flag)
        type(stdlib_state), intent(in) :: err
        integer, intent(in) :: flag
        state_gt_flag = err%ierr > flag
     end function state_gt_flag
     elemental logical function state_ge_flag(err,flag)
        type(stdlib_state), intent(in) :: err
        integer, intent(in) :: flag
        state_ge_flag = err%ierr >= flag
     end function state_ge_flag
     elemental logical function flag_gt_state(flag,err)
        integer, intent(in) :: flag
        type(stdlib_state), intent(in) :: err
        flag_gt_state = err%ierr > flag
     end function flag_gt_state
     elemental logical function flag_ge_state(flag,err)
        integer, intent(in) :: flag
        type(stdlib_state), intent(in) :: err
        flag_ge_state = err%ierr >= flag
     end function flag_ge_state

end module stdlib_linalg_state
