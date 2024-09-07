! Vector norms
module stdlib_linalg_norms
     use stdlib_linalg_constants
     use stdlib_linalg_blas,only:nrm2
     use stdlib_linalg_lapack,only:lange
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private
     
     public :: norm

     character(*),parameter :: this = 'norm'
     
     !> List of internal norm flags
     integer(ilp),parameter :: NORM_ONE = 1
     integer(ilp),parameter :: NORM_TWO = 2
     integer(ilp),parameter :: NORM_P = 3
     integer(ilp),parameter :: NORM_INF = +huge(0_ilp) ! infinity norm
     integer(ilp),parameter :: NORM_MINUSINF = -huge(0_ilp)
     
     !> Vector norm interface
     interface norm
        !> Scalar norms: real(sp)
        module procedure stdlib_linalg_norm_1D_order_s
        module procedure stdlib_linalg_norm_1D_order_err_s
        module procedure stdlib_linalg_norm_2D_order_s
        module procedure stdlib_linalg_norm_2D_order_err_s
        module procedure stdlib_linalg_norm_3D_order_s
        module procedure stdlib_linalg_norm_3D_order_err_s
        module procedure stdlib_linalg_norm_4D_order_s
        module procedure stdlib_linalg_norm_4D_order_err_s
        module procedure stdlib_linalg_norm_5D_order_s
        module procedure stdlib_linalg_norm_5D_order_err_s
        module procedure stdlib_linalg_norm_6D_order_s
        module procedure stdlib_linalg_norm_6D_order_err_s
        module procedure stdlib_linalg_norm_7D_order_s
        module procedure stdlib_linalg_norm_7D_order_err_s
        module procedure stdlib_linalg_norm_8D_order_s
        module procedure stdlib_linalg_norm_8D_order_err_s
        module procedure stdlib_linalg_norm_9D_order_s
        module procedure stdlib_linalg_norm_9D_order_err_s
        module procedure stdlib_linalg_norm_10D_order_s
        module procedure stdlib_linalg_norm_10D_order_err_s
        module procedure stdlib_linalg_norm_11D_order_s
        module procedure stdlib_linalg_norm_11D_order_err_s
        module procedure stdlib_linalg_norm_12D_order_s
        module procedure stdlib_linalg_norm_12D_order_err_s
        module procedure stdlib_linalg_norm_13D_order_s
        module procedure stdlib_linalg_norm_13D_order_err_s
        module procedure stdlib_linalg_norm_14D_order_s
        module procedure stdlib_linalg_norm_14D_order_err_s
        module procedure stdlib_linalg_norm_15D_order_s
        module procedure stdlib_linalg_norm_15D_order_err_s
        !> Array norms: real(sp)
        module procedure stdlib_linalg_norm_2D_to_1D_s
        module procedure stdlib_linalg_norm_2D_to_1D_err_s
        module procedure stdlib_linalg_norm_3D_to_2D_s
        module procedure stdlib_linalg_norm_3D_to_2D_err_s
        module procedure stdlib_linalg_norm_4D_to_3D_s
        module procedure stdlib_linalg_norm_4D_to_3D_err_s
        module procedure stdlib_linalg_norm_5D_to_4D_s
        module procedure stdlib_linalg_norm_5D_to_4D_err_s
        module procedure stdlib_linalg_norm_6D_to_5D_s
        module procedure stdlib_linalg_norm_6D_to_5D_err_s
        module procedure stdlib_linalg_norm_7D_to_6D_s
        module procedure stdlib_linalg_norm_7D_to_6D_err_s
        module procedure stdlib_linalg_norm_8D_to_7D_s
        module procedure stdlib_linalg_norm_8D_to_7D_err_s
        module procedure stdlib_linalg_norm_9D_to_8D_s
        module procedure stdlib_linalg_norm_9D_to_8D_err_s
        module procedure stdlib_linalg_norm_10D_to_9D_s
        module procedure stdlib_linalg_norm_10D_to_9D_err_s
        module procedure stdlib_linalg_norm_11D_to_10D_s
        module procedure stdlib_linalg_norm_11D_to_10D_err_s
        module procedure stdlib_linalg_norm_12D_to_11D_s
        module procedure stdlib_linalg_norm_12D_to_11D_err_s
        module procedure stdlib_linalg_norm_13D_to_12D_s
        module procedure stdlib_linalg_norm_13D_to_12D_err_s
        module procedure stdlib_linalg_norm_14D_to_13D_s
        module procedure stdlib_linalg_norm_14D_to_13D_err_s
        module procedure stdlib_linalg_norm_15D_to_14D_s
        module procedure stdlib_linalg_norm_15D_to_14D_err_s
        !> Scalar norms: real(dp)
        module procedure stdlib_linalg_norm_1D_order_d
        module procedure stdlib_linalg_norm_1D_order_err_d
        module procedure stdlib_linalg_norm_2D_order_d
        module procedure stdlib_linalg_norm_2D_order_err_d
        module procedure stdlib_linalg_norm_3D_order_d
        module procedure stdlib_linalg_norm_3D_order_err_d
        module procedure stdlib_linalg_norm_4D_order_d
        module procedure stdlib_linalg_norm_4D_order_err_d
        module procedure stdlib_linalg_norm_5D_order_d
        module procedure stdlib_linalg_norm_5D_order_err_d
        module procedure stdlib_linalg_norm_6D_order_d
        module procedure stdlib_linalg_norm_6D_order_err_d
        module procedure stdlib_linalg_norm_7D_order_d
        module procedure stdlib_linalg_norm_7D_order_err_d
        module procedure stdlib_linalg_norm_8D_order_d
        module procedure stdlib_linalg_norm_8D_order_err_d
        module procedure stdlib_linalg_norm_9D_order_d
        module procedure stdlib_linalg_norm_9D_order_err_d
        module procedure stdlib_linalg_norm_10D_order_d
        module procedure stdlib_linalg_norm_10D_order_err_d
        module procedure stdlib_linalg_norm_11D_order_d
        module procedure stdlib_linalg_norm_11D_order_err_d
        module procedure stdlib_linalg_norm_12D_order_d
        module procedure stdlib_linalg_norm_12D_order_err_d
        module procedure stdlib_linalg_norm_13D_order_d
        module procedure stdlib_linalg_norm_13D_order_err_d
        module procedure stdlib_linalg_norm_14D_order_d
        module procedure stdlib_linalg_norm_14D_order_err_d
        module procedure stdlib_linalg_norm_15D_order_d
        module procedure stdlib_linalg_norm_15D_order_err_d
        !> Array norms: real(dp)
        module procedure stdlib_linalg_norm_2D_to_1D_d
        module procedure stdlib_linalg_norm_2D_to_1D_err_d
        module procedure stdlib_linalg_norm_3D_to_2D_d
        module procedure stdlib_linalg_norm_3D_to_2D_err_d
        module procedure stdlib_linalg_norm_4D_to_3D_d
        module procedure stdlib_linalg_norm_4D_to_3D_err_d
        module procedure stdlib_linalg_norm_5D_to_4D_d
        module procedure stdlib_linalg_norm_5D_to_4D_err_d
        module procedure stdlib_linalg_norm_6D_to_5D_d
        module procedure stdlib_linalg_norm_6D_to_5D_err_d
        module procedure stdlib_linalg_norm_7D_to_6D_d
        module procedure stdlib_linalg_norm_7D_to_6D_err_d
        module procedure stdlib_linalg_norm_8D_to_7D_d
        module procedure stdlib_linalg_norm_8D_to_7D_err_d
        module procedure stdlib_linalg_norm_9D_to_8D_d
        module procedure stdlib_linalg_norm_9D_to_8D_err_d
        module procedure stdlib_linalg_norm_10D_to_9D_d
        module procedure stdlib_linalg_norm_10D_to_9D_err_d
        module procedure stdlib_linalg_norm_11D_to_10D_d
        module procedure stdlib_linalg_norm_11D_to_10D_err_d
        module procedure stdlib_linalg_norm_12D_to_11D_d
        module procedure stdlib_linalg_norm_12D_to_11D_err_d
        module procedure stdlib_linalg_norm_13D_to_12D_d
        module procedure stdlib_linalg_norm_13D_to_12D_err_d
        module procedure stdlib_linalg_norm_14D_to_13D_d
        module procedure stdlib_linalg_norm_14D_to_13D_err_d
        module procedure stdlib_linalg_norm_15D_to_14D_d
        module procedure stdlib_linalg_norm_15D_to_14D_err_d
        !> Scalar norms: real(qp)
        module procedure stdlib_linalg_norm_1D_order_q
        module procedure stdlib_linalg_norm_1D_order_err_q
        module procedure stdlib_linalg_norm_2D_order_q
        module procedure stdlib_linalg_norm_2D_order_err_q
        module procedure stdlib_linalg_norm_3D_order_q
        module procedure stdlib_linalg_norm_3D_order_err_q
        module procedure stdlib_linalg_norm_4D_order_q
        module procedure stdlib_linalg_norm_4D_order_err_q
        module procedure stdlib_linalg_norm_5D_order_q
        module procedure stdlib_linalg_norm_5D_order_err_q
        module procedure stdlib_linalg_norm_6D_order_q
        module procedure stdlib_linalg_norm_6D_order_err_q
        module procedure stdlib_linalg_norm_7D_order_q
        module procedure stdlib_linalg_norm_7D_order_err_q
        module procedure stdlib_linalg_norm_8D_order_q
        module procedure stdlib_linalg_norm_8D_order_err_q
        module procedure stdlib_linalg_norm_9D_order_q
        module procedure stdlib_linalg_norm_9D_order_err_q
        module procedure stdlib_linalg_norm_10D_order_q
        module procedure stdlib_linalg_norm_10D_order_err_q
        module procedure stdlib_linalg_norm_11D_order_q
        module procedure stdlib_linalg_norm_11D_order_err_q
        module procedure stdlib_linalg_norm_12D_order_q
        module procedure stdlib_linalg_norm_12D_order_err_q
        module procedure stdlib_linalg_norm_13D_order_q
        module procedure stdlib_linalg_norm_13D_order_err_q
        module procedure stdlib_linalg_norm_14D_order_q
        module procedure stdlib_linalg_norm_14D_order_err_q
        module procedure stdlib_linalg_norm_15D_order_q
        module procedure stdlib_linalg_norm_15D_order_err_q
        !> Array norms: real(qp)
        module procedure stdlib_linalg_norm_2D_to_1D_q
        module procedure stdlib_linalg_norm_2D_to_1D_err_q
        module procedure stdlib_linalg_norm_3D_to_2D_q
        module procedure stdlib_linalg_norm_3D_to_2D_err_q
        module procedure stdlib_linalg_norm_4D_to_3D_q
        module procedure stdlib_linalg_norm_4D_to_3D_err_q
        module procedure stdlib_linalg_norm_5D_to_4D_q
        module procedure stdlib_linalg_norm_5D_to_4D_err_q
        module procedure stdlib_linalg_norm_6D_to_5D_q
        module procedure stdlib_linalg_norm_6D_to_5D_err_q
        module procedure stdlib_linalg_norm_7D_to_6D_q
        module procedure stdlib_linalg_norm_7D_to_6D_err_q
        module procedure stdlib_linalg_norm_8D_to_7D_q
        module procedure stdlib_linalg_norm_8D_to_7D_err_q
        module procedure stdlib_linalg_norm_9D_to_8D_q
        module procedure stdlib_linalg_norm_9D_to_8D_err_q
        module procedure stdlib_linalg_norm_10D_to_9D_q
        module procedure stdlib_linalg_norm_10D_to_9D_err_q
        module procedure stdlib_linalg_norm_11D_to_10D_q
        module procedure stdlib_linalg_norm_11D_to_10D_err_q
        module procedure stdlib_linalg_norm_12D_to_11D_q
        module procedure stdlib_linalg_norm_12D_to_11D_err_q
        module procedure stdlib_linalg_norm_13D_to_12D_q
        module procedure stdlib_linalg_norm_13D_to_12D_err_q
        module procedure stdlib_linalg_norm_14D_to_13D_q
        module procedure stdlib_linalg_norm_14D_to_13D_err_q
        module procedure stdlib_linalg_norm_15D_to_14D_q
        module procedure stdlib_linalg_norm_15D_to_14D_err_q
        !> Scalar norms: complex(sp)
        module procedure stdlib_linalg_norm_1D_order_c
        module procedure stdlib_linalg_norm_1D_order_err_c
        module procedure stdlib_linalg_norm_2D_order_c
        module procedure stdlib_linalg_norm_2D_order_err_c
        module procedure stdlib_linalg_norm_3D_order_c
        module procedure stdlib_linalg_norm_3D_order_err_c
        module procedure stdlib_linalg_norm_4D_order_c
        module procedure stdlib_linalg_norm_4D_order_err_c
        module procedure stdlib_linalg_norm_5D_order_c
        module procedure stdlib_linalg_norm_5D_order_err_c
        module procedure stdlib_linalg_norm_6D_order_c
        module procedure stdlib_linalg_norm_6D_order_err_c
        module procedure stdlib_linalg_norm_7D_order_c
        module procedure stdlib_linalg_norm_7D_order_err_c
        module procedure stdlib_linalg_norm_8D_order_c
        module procedure stdlib_linalg_norm_8D_order_err_c
        module procedure stdlib_linalg_norm_9D_order_c
        module procedure stdlib_linalg_norm_9D_order_err_c
        module procedure stdlib_linalg_norm_10D_order_c
        module procedure stdlib_linalg_norm_10D_order_err_c
        module procedure stdlib_linalg_norm_11D_order_c
        module procedure stdlib_linalg_norm_11D_order_err_c
        module procedure stdlib_linalg_norm_12D_order_c
        module procedure stdlib_linalg_norm_12D_order_err_c
        module procedure stdlib_linalg_norm_13D_order_c
        module procedure stdlib_linalg_norm_13D_order_err_c
        module procedure stdlib_linalg_norm_14D_order_c
        module procedure stdlib_linalg_norm_14D_order_err_c
        module procedure stdlib_linalg_norm_15D_order_c
        module procedure stdlib_linalg_norm_15D_order_err_c
        !> Array norms: complex(sp)
        module procedure stdlib_linalg_norm_2D_to_1D_c
        module procedure stdlib_linalg_norm_2D_to_1D_err_c
        module procedure stdlib_linalg_norm_3D_to_2D_c
        module procedure stdlib_linalg_norm_3D_to_2D_err_c
        module procedure stdlib_linalg_norm_4D_to_3D_c
        module procedure stdlib_linalg_norm_4D_to_3D_err_c
        module procedure stdlib_linalg_norm_5D_to_4D_c
        module procedure stdlib_linalg_norm_5D_to_4D_err_c
        module procedure stdlib_linalg_norm_6D_to_5D_c
        module procedure stdlib_linalg_norm_6D_to_5D_err_c
        module procedure stdlib_linalg_norm_7D_to_6D_c
        module procedure stdlib_linalg_norm_7D_to_6D_err_c
        module procedure stdlib_linalg_norm_8D_to_7D_c
        module procedure stdlib_linalg_norm_8D_to_7D_err_c
        module procedure stdlib_linalg_norm_9D_to_8D_c
        module procedure stdlib_linalg_norm_9D_to_8D_err_c
        module procedure stdlib_linalg_norm_10D_to_9D_c
        module procedure stdlib_linalg_norm_10D_to_9D_err_c
        module procedure stdlib_linalg_norm_11D_to_10D_c
        module procedure stdlib_linalg_norm_11D_to_10D_err_c
        module procedure stdlib_linalg_norm_12D_to_11D_c
        module procedure stdlib_linalg_norm_12D_to_11D_err_c
        module procedure stdlib_linalg_norm_13D_to_12D_c
        module procedure stdlib_linalg_norm_13D_to_12D_err_c
        module procedure stdlib_linalg_norm_14D_to_13D_c
        module procedure stdlib_linalg_norm_14D_to_13D_err_c
        module procedure stdlib_linalg_norm_15D_to_14D_c
        module procedure stdlib_linalg_norm_15D_to_14D_err_c
        !> Scalar norms: complex(dp)
        module procedure stdlib_linalg_norm_1D_order_z
        module procedure stdlib_linalg_norm_1D_order_err_z
        module procedure stdlib_linalg_norm_2D_order_z
        module procedure stdlib_linalg_norm_2D_order_err_z
        module procedure stdlib_linalg_norm_3D_order_z
        module procedure stdlib_linalg_norm_3D_order_err_z
        module procedure stdlib_linalg_norm_4D_order_z
        module procedure stdlib_linalg_norm_4D_order_err_z
        module procedure stdlib_linalg_norm_5D_order_z
        module procedure stdlib_linalg_norm_5D_order_err_z
        module procedure stdlib_linalg_norm_6D_order_z
        module procedure stdlib_linalg_norm_6D_order_err_z
        module procedure stdlib_linalg_norm_7D_order_z
        module procedure stdlib_linalg_norm_7D_order_err_z
        module procedure stdlib_linalg_norm_8D_order_z
        module procedure stdlib_linalg_norm_8D_order_err_z
        module procedure stdlib_linalg_norm_9D_order_z
        module procedure stdlib_linalg_norm_9D_order_err_z
        module procedure stdlib_linalg_norm_10D_order_z
        module procedure stdlib_linalg_norm_10D_order_err_z
        module procedure stdlib_linalg_norm_11D_order_z
        module procedure stdlib_linalg_norm_11D_order_err_z
        module procedure stdlib_linalg_norm_12D_order_z
        module procedure stdlib_linalg_norm_12D_order_err_z
        module procedure stdlib_linalg_norm_13D_order_z
        module procedure stdlib_linalg_norm_13D_order_err_z
        module procedure stdlib_linalg_norm_14D_order_z
        module procedure stdlib_linalg_norm_14D_order_err_z
        module procedure stdlib_linalg_norm_15D_order_z
        module procedure stdlib_linalg_norm_15D_order_err_z
        !> Array norms: complex(dp)
        module procedure stdlib_linalg_norm_2D_to_1D_z
        module procedure stdlib_linalg_norm_2D_to_1D_err_z
        module procedure stdlib_linalg_norm_3D_to_2D_z
        module procedure stdlib_linalg_norm_3D_to_2D_err_z
        module procedure stdlib_linalg_norm_4D_to_3D_z
        module procedure stdlib_linalg_norm_4D_to_3D_err_z
        module procedure stdlib_linalg_norm_5D_to_4D_z
        module procedure stdlib_linalg_norm_5D_to_4D_err_z
        module procedure stdlib_linalg_norm_6D_to_5D_z
        module procedure stdlib_linalg_norm_6D_to_5D_err_z
        module procedure stdlib_linalg_norm_7D_to_6D_z
        module procedure stdlib_linalg_norm_7D_to_6D_err_z
        module procedure stdlib_linalg_norm_8D_to_7D_z
        module procedure stdlib_linalg_norm_8D_to_7D_err_z
        module procedure stdlib_linalg_norm_9D_to_8D_z
        module procedure stdlib_linalg_norm_9D_to_8D_err_z
        module procedure stdlib_linalg_norm_10D_to_9D_z
        module procedure stdlib_linalg_norm_10D_to_9D_err_z
        module procedure stdlib_linalg_norm_11D_to_10D_z
        module procedure stdlib_linalg_norm_11D_to_10D_err_z
        module procedure stdlib_linalg_norm_12D_to_11D_z
        module procedure stdlib_linalg_norm_12D_to_11D_err_z
        module procedure stdlib_linalg_norm_13D_to_12D_z
        module procedure stdlib_linalg_norm_13D_to_12D_err_z
        module procedure stdlib_linalg_norm_14D_to_13D_z
        module procedure stdlib_linalg_norm_14D_to_13D_err_z
        module procedure stdlib_linalg_norm_15D_to_14D_z
        module procedure stdlib_linalg_norm_15D_to_14D_err_z
        !> Scalar norms: complex(qp)
        module procedure stdlib_linalg_norm_1D_order_w
        module procedure stdlib_linalg_norm_1D_order_err_w
        module procedure stdlib_linalg_norm_2D_order_w
        module procedure stdlib_linalg_norm_2D_order_err_w
        module procedure stdlib_linalg_norm_3D_order_w
        module procedure stdlib_linalg_norm_3D_order_err_w
        module procedure stdlib_linalg_norm_4D_order_w
        module procedure stdlib_linalg_norm_4D_order_err_w
        module procedure stdlib_linalg_norm_5D_order_w
        module procedure stdlib_linalg_norm_5D_order_err_w
        module procedure stdlib_linalg_norm_6D_order_w
        module procedure stdlib_linalg_norm_6D_order_err_w
        module procedure stdlib_linalg_norm_7D_order_w
        module procedure stdlib_linalg_norm_7D_order_err_w
        module procedure stdlib_linalg_norm_8D_order_w
        module procedure stdlib_linalg_norm_8D_order_err_w
        module procedure stdlib_linalg_norm_9D_order_w
        module procedure stdlib_linalg_norm_9D_order_err_w
        module procedure stdlib_linalg_norm_10D_order_w
        module procedure stdlib_linalg_norm_10D_order_err_w
        module procedure stdlib_linalg_norm_11D_order_w
        module procedure stdlib_linalg_norm_11D_order_err_w
        module procedure stdlib_linalg_norm_12D_order_w
        module procedure stdlib_linalg_norm_12D_order_err_w
        module procedure stdlib_linalg_norm_13D_order_w
        module procedure stdlib_linalg_norm_13D_order_err_w
        module procedure stdlib_linalg_norm_14D_order_w
        module procedure stdlib_linalg_norm_14D_order_err_w
        module procedure stdlib_linalg_norm_15D_order_w
        module procedure stdlib_linalg_norm_15D_order_err_w
        !> Array norms: complex(qp)
        module procedure stdlib_linalg_norm_2D_to_1D_w
        module procedure stdlib_linalg_norm_2D_to_1D_err_w
        module procedure stdlib_linalg_norm_3D_to_2D_w
        module procedure stdlib_linalg_norm_3D_to_2D_err_w
        module procedure stdlib_linalg_norm_4D_to_3D_w
        module procedure stdlib_linalg_norm_4D_to_3D_err_w
        module procedure stdlib_linalg_norm_5D_to_4D_w
        module procedure stdlib_linalg_norm_5D_to_4D_err_w
        module procedure stdlib_linalg_norm_6D_to_5D_w
        module procedure stdlib_linalg_norm_6D_to_5D_err_w
        module procedure stdlib_linalg_norm_7D_to_6D_w
        module procedure stdlib_linalg_norm_7D_to_6D_err_w
        module procedure stdlib_linalg_norm_8D_to_7D_w
        module procedure stdlib_linalg_norm_8D_to_7D_err_w
        module procedure stdlib_linalg_norm_9D_to_8D_w
        module procedure stdlib_linalg_norm_9D_to_8D_err_w
        module procedure stdlib_linalg_norm_10D_to_9D_w
        module procedure stdlib_linalg_norm_10D_to_9D_err_w
        module procedure stdlib_linalg_norm_11D_to_10D_w
        module procedure stdlib_linalg_norm_11D_to_10D_err_w
        module procedure stdlib_linalg_norm_12D_to_11D_w
        module procedure stdlib_linalg_norm_12D_to_11D_err_w
        module procedure stdlib_linalg_norm_13D_to_12D_w
        module procedure stdlib_linalg_norm_13D_to_12D_err_w
        module procedure stdlib_linalg_norm_14D_to_13D_w
        module procedure stdlib_linalg_norm_14D_to_13D_err_w
        module procedure stdlib_linalg_norm_15D_to_14D_w
        module procedure stdlib_linalg_norm_15D_to_14D_err_w
     end interface norm
     
     interface parse_norm_type
        module procedure parse_norm_type_integer
        module procedure parse_norm_type_character
     end interface parse_norm_type
     
     contains
     
     !> Parse norm type from an integer user input
     pure subroutine parse_norm_type_integer(order,norm_type,err)
        !> User input value
        integer(ilp),intent(in) :: order
        !> Return value: norm type
        integer(ilp),intent(out) :: norm_type
        !> State return flag
        type(linalg_state),intent(out) :: err
        
        select case (order)
           case (1_ilp)
               norm_type = NORM_ONE
           case (2_ilp)
               norm_type = NORM_TWO
           case (3_ilp:huge(0_ilp) - 1_ilp)
               norm_type = NORM_P
           case (huge(0_ilp):)
               norm_type = NORM_INF
           case (:-huge(0_ilp))
               norm_type = NORM_MINUSINF
           
           case default
               norm_type = NORM_ONE
               err = linalg_state(this,LINALG_ERROR,'Input norm type ',order,' is not recognized.')
        end select
        
     end subroutine parse_norm_type_integer

     pure subroutine parse_norm_type_character(order,norm_type,err)
        !> User input value
        character(*),intent(in) :: order
        !> Return value: norm type
        integer(ilp),intent(out) :: norm_type
        !> State return flag
        type(linalg_state),intent(out) :: err
        
        integer(ilp) :: int_order,read_err
        
        select case (order)
           case ('inf','Inf','INF')
              norm_type = NORM_INF
           case ('-inf','-Inf','-INF')
              norm_type = NORM_MINUSINF
           case ('Euclidean','euclidean','EUCLIDEAN')
              norm_type = NORM_TWO
           case default
            
              ! Check if this input can be read as an integer
              read (order,*,iostat=read_err) int_order
              if (read_err /= 0) then
                 ! Cannot read as an integer
                 norm_type = NORM_ONE
                 err = linalg_state(this,LINALG_ERROR,'Input norm type ',order,' is not recognized.')
              else
                 call parse_norm_type_integer(int_order,norm_type,err)
              end if

        end select
        
     end subroutine parse_norm_type_character

    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_1D_order_s(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_1D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_1D_order_err_s(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_1D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_1D_s(a,order,err,nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_1D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_2D_order_s(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_2D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_2D_order_err_s(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_2D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_2D_s(a,order,err,nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_3D_order_s(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_3D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_3D_order_err_s(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_3D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_3D_s(a,order,err,nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_4D_order_s(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_4D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_4D_order_err_s(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_4D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_4D_s(a,order,err,nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_5D_order_s(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_5D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_5D_order_err_s(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_5D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_5D_s(a,order,err,nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_6D_order_s(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_6D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_6D_order_err_s(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_6D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_6D_s(a,order,err,nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_7D_order_s(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_7D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_7D_order_err_s(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_7D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_7D_s(a,order,err,nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_8D_order_s(a,order) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_8D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_8D_order_err_s(a,order,err) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_8D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_8D_s(a,order,err,nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_9D_order_s(a,order) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_9D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_9D_order_err_s(a,order,err) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_9D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_9D_s(a,order,err,nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_10D_order_s(a,order) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_10D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_10D_order_err_s(a,order,err) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_10D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_10D_s(a,order,err,nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_11D_order_s(a,order) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_11D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_11D_order_err_s(a,order,err) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_11D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_11D_s(a,order,err,nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_12D_order_s(a,order) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_12D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_12D_order_err_s(a,order,err) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_12D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_12D_s(a,order,err,nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_13D_order_s(a,order) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_13D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_13D_order_err_s(a,order,err) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_13D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_13D_s(a,order,err,nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_14D_order_s(a,order) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_14D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_14D_order_err_s(a,order,err) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_14D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_14D_s(a,order,err,nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_15D_order_s(a,order) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_15D_s(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_s
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_15D_order_err_s(a,order,err) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_15D_s(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_err_s
    
    ! Internal implementation
    pure subroutine norm_15D_s(a,order,err,nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_1D_order_d(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_1D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_1D_order_err_d(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_1D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_1D_d(a,order,err,nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_1D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_2D_order_d(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_2D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_2D_order_err_d(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_2D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_2D_d(a,order,err,nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_3D_order_d(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_3D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_3D_order_err_d(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_3D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_3D_d(a,order,err,nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_4D_order_d(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_4D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_4D_order_err_d(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_4D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_4D_d(a,order,err,nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_5D_order_d(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_5D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_5D_order_err_d(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_5D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_5D_d(a,order,err,nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_6D_order_d(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_6D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_6D_order_err_d(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_6D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_6D_d(a,order,err,nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_7D_order_d(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_7D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_7D_order_err_d(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_7D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_7D_d(a,order,err,nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_8D_order_d(a,order) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_8D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_8D_order_err_d(a,order,err) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_8D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_8D_d(a,order,err,nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_9D_order_d(a,order) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_9D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_9D_order_err_d(a,order,err) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_9D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_9D_d(a,order,err,nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_10D_order_d(a,order) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_10D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_10D_order_err_d(a,order,err) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_10D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_10D_d(a,order,err,nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_11D_order_d(a,order) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_11D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_11D_order_err_d(a,order,err) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_11D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_11D_d(a,order,err,nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_12D_order_d(a,order) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_12D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_12D_order_err_d(a,order,err) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_12D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_12D_d(a,order,err,nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_13D_order_d(a,order) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_13D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_13D_order_err_d(a,order,err) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_13D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_13D_d(a,order,err,nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_14D_order_d(a,order) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_14D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_14D_order_err_d(a,order,err) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_14D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_14D_d(a,order,err,nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_15D_order_d(a,order) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_15D_d(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_d
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_15D_order_err_d(a,order,err) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_15D_d(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_err_d
    
    ! Internal implementation
    pure subroutine norm_15D_d(a,order,err,nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_1D_order_q(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_1D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_1D_order_err_q(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_1D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_1D_q(a,order,err,nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_1D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_2D_order_q(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_2D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_2D_order_err_q(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_2D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_2D_q(a,order,err,nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_3D_order_q(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_3D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_3D_order_err_q(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_3D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_3D_q(a,order,err,nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_4D_order_q(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_4D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_4D_order_err_q(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_4D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_4D_q(a,order,err,nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_5D_order_q(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_5D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_5D_order_err_q(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_5D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_5D_q(a,order,err,nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_6D_order_q(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_6D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_6D_order_err_q(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_6D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_6D_q(a,order,err,nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_7D_order_q(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_7D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_7D_order_err_q(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_7D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_7D_q(a,order,err,nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_8D_order_q(a,order) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_8D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_8D_order_err_q(a,order,err) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_8D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_8D_q(a,order,err,nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_9D_order_q(a,order) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_9D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_9D_order_err_q(a,order,err) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_9D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_9D_q(a,order,err,nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_10D_order_q(a,order) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_10D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_10D_order_err_q(a,order,err) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_10D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_10D_q(a,order,err,nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_11D_order_q(a,order) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_11D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_11D_order_err_q(a,order,err) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_11D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_11D_q(a,order,err,nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_12D_order_q(a,order) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_12D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_12D_order_err_q(a,order,err) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_12D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_12D_q(a,order,err,nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_13D_order_q(a,order) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_13D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_13D_order_err_q(a,order,err) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_13D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_13D_q(a,order,err,nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_14D_order_q(a,order) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_14D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_14D_order_err_q(a,order,err) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_14D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_14D_q(a,order,err,nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_15D_order_q(a,order) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_15D_q(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_q
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_15D_order_err_q(a,order,err) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_15D_q(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_err_q
    
    ! Internal implementation
    pure subroutine norm_15D_q(a,order,err,nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a**2))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_1D_order_c(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_1D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_1D_order_err_c(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_1D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_1D_c(a,order,err,nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_1D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_2D_order_c(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_2D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_2D_order_err_c(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_2D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_2D_c(a,order,err,nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_3D_order_c(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_3D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_3D_order_err_c(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_3D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_3D_c(a,order,err,nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_4D_order_c(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_4D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_4D_order_err_c(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_4D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_4D_c(a,order,err,nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_5D_order_c(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_5D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_5D_order_err_c(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_5D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_5D_c(a,order,err,nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_6D_order_c(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_6D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_6D_order_err_c(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_6D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_6D_c(a,order,err,nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_7D_order_c(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_7D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_7D_order_err_c(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_7D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_7D_c(a,order,err,nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_8D_order_c(a,order) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_8D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_8D_order_err_c(a,order,err) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_8D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_8D_c(a,order,err,nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_9D_order_c(a,order) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_9D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_9D_order_err_c(a,order,err) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_9D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_9D_c(a,order,err,nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_10D_order_c(a,order) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_10D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_10D_order_err_c(a,order,err) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_10D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_10D_c(a,order,err,nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_11D_order_c(a,order) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_11D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_11D_order_err_c(a,order,err) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_11D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_11D_c(a,order,err,nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_12D_order_c(a,order) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_12D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_12D_order_err_c(a,order,err) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_12D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_12D_c(a,order,err,nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_13D_order_c(a,order) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_13D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_13D_order_err_c(a,order,err) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_13D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_13D_c(a,order,err,nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_14D_order_c(a,order) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_14D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_14D_order_err_c(a,order,err) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_14D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_14D_c(a,order,err,nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_15D_order_c(a,order) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_15D_c(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_c
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_15D_order_err_c(a,order,err) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm
                
        call norm_15D_c(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_err_c
    
    ! Internal implementation
    pure subroutine norm_15D_c(a,order,err,nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_1D_order_z(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_1D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_1D_order_err_z(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_1D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_1D_z(a,order,err,nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_1D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_2D_order_z(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_2D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_2D_order_err_z(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_2D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_2D_z(a,order,err,nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_3D_order_z(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_3D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_3D_order_err_z(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_3D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_3D_z(a,order,err,nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_4D_order_z(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_4D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_4D_order_err_z(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_4D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_4D_z(a,order,err,nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_5D_order_z(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_5D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_5D_order_err_z(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_5D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_5D_z(a,order,err,nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_6D_order_z(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_6D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_6D_order_err_z(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_6D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_6D_z(a,order,err,nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_7D_order_z(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_7D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_7D_order_err_z(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_7D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_7D_z(a,order,err,nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_8D_order_z(a,order) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_8D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_8D_order_err_z(a,order,err) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_8D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_8D_z(a,order,err,nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_9D_order_z(a,order) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_9D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_9D_order_err_z(a,order,err) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_9D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_9D_z(a,order,err,nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_10D_order_z(a,order) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_10D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_10D_order_err_z(a,order,err) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_10D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_10D_z(a,order,err,nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_11D_order_z(a,order) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_11D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_11D_order_err_z(a,order,err) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_11D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_11D_z(a,order,err,nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_12D_order_z(a,order) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_12D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_12D_order_err_z(a,order,err) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_12D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_12D_z(a,order,err,nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_13D_order_z(a,order) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_13D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_13D_order_err_z(a,order,err) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_13D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_13D_z(a,order,err,nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_14D_order_z(a,order) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_14D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_14D_order_err_z(a,order,err) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_14D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_14D_z(a,order,err,nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_15D_order_z(a,order) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_15D_z(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_z
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_15D_order_err_z(a,order,err) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm
                
        call norm_15D_z(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_err_z
    
    ! Internal implementation
    pure subroutine norm_15D_z(a,order,err,nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_1D_order_w(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_1D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_1D_order_err_w(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_1D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_1D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_1D_w(a,order,err,nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_1D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_2D_order_w(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_2D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_2D_order_err_w(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_2D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_2D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_2D_w(a,order,err,nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_3D_order_w(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_3D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_3D_order_err_w(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_3D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_3D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_3D_w(a,order,err,nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_4D_order_w(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_4D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_4D_order_err_w(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_4D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_4D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_4D_w(a,order,err,nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_5D_order_w(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_5D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_5D_order_err_w(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_5D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_5D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_5D_w(a,order,err,nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_6D_order_w(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_6D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_6D_order_err_w(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_6D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_6D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_6D_w(a,order,err,nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_7D_order_w(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_7D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_7D_order_err_w(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_7D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_7D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_7D_w(a,order,err,nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_8D_order_w(a,order) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_8D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_8D_order_err_w(a,order,err) result(nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_8D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_8D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_8D_w(a,order,err,nrm)
        !> Input 8-d matrix a(:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_9D_order_w(a,order) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_9D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_9D_order_err_w(a,order,err) result(nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_9D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_9D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_9D_w(a,order,err,nrm)
        !> Input 9-d matrix a(:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_10D_order_w(a,order) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_10D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_10D_order_err_w(a,order,err) result(nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_10D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_10D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_10D_w(a,order,err,nrm)
        !> Input 10-d matrix a(:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_11D_order_w(a,order) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_11D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_11D_order_err_w(a,order,err) result(nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_11D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_11D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_11D_w(a,order,err,nrm)
        !> Input 11-d matrix a(:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_12D_order_w(a,order) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_12D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_12D_order_err_w(a,order,err) result(nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_12D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_12D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_12D_w(a,order,err,nrm)
        !> Input 12-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_13D_order_w(a,order) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_13D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_13D_order_err_w(a,order,err) result(nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_13D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_13D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_13D_w(a,order,err,nrm)
        !> Input 13-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_14D_order_w(a,order) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_14D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_14D_order_err_w(a,order,err) result(nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_14D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_14D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_14D_w(a,order,err,nrm)
        !> Input 14-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function stdlib_linalg_norm_15D_order_w(a,order) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_15D_w(a,order=order,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_w
    
    ! Function interface with output error
    ! Pure function interface, with order specification
    function stdlib_linalg_norm_15D_order_err_w(a,order,err) result(nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm
                
        call norm_15D_w(a,order=order,err=err,nrm=nrm)
        
    end function stdlib_linalg_norm_15D_order_err_w
    
    ! Internal implementation
    pure subroutine norm_15D_w(a,order,err,nrm)
        !> Input 15-d matrix a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        
        type(linalg_state) :: err_
        
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a)))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_w

    !===============================================
    ! Norms : any rank to rank-1, with DIM specifier
    !===============================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_2D_to_1D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_2D_to_1D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_err_s
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 2) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_to_1D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_3D_to_2D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_3D_to_2D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_err_s
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 3) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_to_2D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_4D_to_3D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_4D_to_3D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_err_s
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 4) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_to_3D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_5D_to_4D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_5D_to_4D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_err_s
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 5) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_to_4D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_6D_to_5D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_6D_to_5D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_err_s
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 6) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_to_5D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_7D_to_6D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_7D_to_6D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_err_s
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 7) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_to_6D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_8D_to_7D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_8D_to_7D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_err_s
    
    ! Internal implementation
    pure subroutine norm_8D_to_7D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 8) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_to_7D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_9D_to_8D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_9D_to_8D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_err_s
    
    ! Internal implementation
    pure subroutine norm_9D_to_8D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 9) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_to_8D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_10D_to_9D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_10D_to_9D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_err_s
    
    ! Internal implementation
    pure subroutine norm_10D_to_9D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 10) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_to_9D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_11D_to_10D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_11D_to_10D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_err_s
    
    ! Internal implementation
    pure subroutine norm_11D_to_10D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))  &
            &   
type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 11) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_to_10D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_12D_to_11D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_12D_to_11D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_err_s
    
    ! Internal implementation
    pure subroutine norm_12D_to_11D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 12) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_to_11D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_13D_to_12D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_13D_to_12D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_err_s
    
    ! Internal implementation
    pure subroutine norm_13D_to_12D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 13) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_to_12D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_14D_to_13D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_14D_to_13D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_err_s
    
    ! Internal implementation
    pure subroutine norm_14D_to_13D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 14) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_to_13D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_15D_to_14D_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_s(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_s

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_15D_to_14D_err_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_s(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_err_s
    
    ! Internal implementation
    pure subroutine norm_15D_to_14D_s(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 15) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_to_14D_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_2D_to_1D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_2D_to_1D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_err_d
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 2) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_to_1D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_3D_to_2D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_3D_to_2D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_err_d
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 3) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_to_2D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_4D_to_3D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_4D_to_3D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_err_d
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 4) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_to_3D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_5D_to_4D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_5D_to_4D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_err_d
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 5) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_to_4D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_6D_to_5D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_6D_to_5D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_err_d
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 6) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_to_5D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_7D_to_6D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_7D_to_6D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_err_d
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 7) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_to_6D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_8D_to_7D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_8D_to_7D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_err_d
    
    ! Internal implementation
    pure subroutine norm_8D_to_7D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 8) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_to_7D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_9D_to_8D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_9D_to_8D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_err_d
    
    ! Internal implementation
    pure subroutine norm_9D_to_8D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 9) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_to_8D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_10D_to_9D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_10D_to_9D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_err_d
    
    ! Internal implementation
    pure subroutine norm_10D_to_9D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 10) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_to_9D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_11D_to_10D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_11D_to_10D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_err_d
    
    ! Internal implementation
    pure subroutine norm_11D_to_10D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))  &
            &   
type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 11) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_to_10D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_12D_to_11D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_12D_to_11D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_err_d
    
    ! Internal implementation
    pure subroutine norm_12D_to_11D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 12) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_to_11D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_13D_to_12D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_13D_to_12D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_err_d
    
    ! Internal implementation
    pure subroutine norm_13D_to_12D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 13) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_to_12D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_14D_to_13D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_14D_to_13D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_err_d
    
    ! Internal implementation
    pure subroutine norm_14D_to_13D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 14) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_to_13D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_15D_to_14D_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_d(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_d

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_15D_to_14D_err_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_d(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_err_d
    
    ! Internal implementation
    pure subroutine norm_15D_to_14D_d(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 15) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_to_14D_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_2D_to_1D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_2D_to_1D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_err_q
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 2) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_to_1D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_3D_to_2D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_3D_to_2D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_err_q
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 3) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_to_2D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_4D_to_3D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_4D_to_3D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_err_q
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 4) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_to_3D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_5D_to_4D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_5D_to_4D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_err_q
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 5) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_to_4D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_6D_to_5D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_6D_to_5D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_err_q
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 6) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_to_5D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_7D_to_6D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_7D_to_6D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_err_q
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 7) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_to_6D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_8D_to_7D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_8D_to_7D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_err_q
    
    ! Internal implementation
    pure subroutine norm_8D_to_7D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 8) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_to_7D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_9D_to_8D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_9D_to_8D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_err_q
    
    ! Internal implementation
    pure subroutine norm_9D_to_8D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 9) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_to_8D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_10D_to_9D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_10D_to_9D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_err_q
    
    ! Internal implementation
    pure subroutine norm_10D_to_9D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 10) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_to_9D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_11D_to_10D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_11D_to_10D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_err_q
    
    ! Internal implementation
    pure subroutine norm_11D_to_10D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))  &
            &   
type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 11) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_to_10D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_12D_to_11D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_12D_to_11D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_err_q
    
    ! Internal implementation
    pure subroutine norm_12D_to_11D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 12) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_to_11D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_13D_to_12D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_13D_to_12D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_err_q
    
    ! Internal implementation
    pure subroutine norm_13D_to_12D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 13) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_to_12D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_14D_to_13D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_14D_to_13D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_err_q
    
    ! Internal implementation
    pure subroutine norm_14D_to_13D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 14) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_to_13D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_15D_to_14D_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_q(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_q

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_15D_to_14D_err_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_q(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_err_q
    
    ! Internal implementation
    pure subroutine norm_15D_to_14D_q(a,order,dim,err,nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 15) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_to_14D_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_2D_to_1D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_2D_to_1D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_err_c
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 2) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_to_1D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_3D_to_2D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_3D_to_2D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_err_c
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 3) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_to_2D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_4D_to_3D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_4D_to_3D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_err_c
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 4) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_to_3D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_5D_to_4D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_5D_to_4D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_err_c
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 5) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_to_4D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_6D_to_5D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_6D_to_5D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_err_c
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 6) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_to_5D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_7D_to_6D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_7D_to_6D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_err_c
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 7) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_to_6D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_8D_to_7D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_8D_to_7D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_err_c
    
    ! Internal implementation
    pure subroutine norm_8D_to_7D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 8) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_to_7D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_9D_to_8D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_9D_to_8D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_err_c
    
    ! Internal implementation
    pure subroutine norm_9D_to_8D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 9) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_to_8D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_10D_to_9D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_10D_to_9D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_err_c
    
    ! Internal implementation
    pure subroutine norm_10D_to_9D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 10) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_to_9D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_11D_to_10D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_11D_to_10D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_err_c
    
    ! Internal implementation
    pure subroutine norm_11D_to_10D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))  &
            &   
type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 11) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_to_10D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_12D_to_11D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_12D_to_11D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_err_c
    
    ! Internal implementation
    pure subroutine norm_12D_to_11D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 12) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_to_11D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_13D_to_12D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_13D_to_12D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_err_c
    
    ! Internal implementation
    pure subroutine norm_13D_to_12D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 13) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_to_12D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_14D_to_13D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_14D_to_13D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_err_c
    
    ! Internal implementation
    pure subroutine norm_14D_to_13D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 14) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_to_13D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_15D_to_14D_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_c(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_c

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_15D_to_14D_err_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_c(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_err_c
    
    ! Internal implementation
    pure subroutine norm_15D_to_14D_c(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 15) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_to_14D_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_2D_to_1D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_2D_to_1D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_err_z
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 2) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_to_1D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_3D_to_2D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_3D_to_2D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_err_z
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 3) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_to_2D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_4D_to_3D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_4D_to_3D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_err_z
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 4) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_to_3D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_5D_to_4D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_5D_to_4D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_err_z
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 5) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_to_4D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_6D_to_5D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_6D_to_5D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_err_z
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 6) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_to_5D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_7D_to_6D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_7D_to_6D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_err_z
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 7) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_to_6D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_8D_to_7D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_8D_to_7D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_err_z
    
    ! Internal implementation
    pure subroutine norm_8D_to_7D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 8) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_to_7D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_9D_to_8D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_9D_to_8D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_err_z
    
    ! Internal implementation
    pure subroutine norm_9D_to_8D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 9) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_to_8D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_10D_to_9D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_10D_to_9D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_err_z
    
    ! Internal implementation
    pure subroutine norm_10D_to_9D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 10) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_to_9D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_11D_to_10D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_11D_to_10D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_err_z
    
    ! Internal implementation
    pure subroutine norm_11D_to_10D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))  &
            &   
type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 11) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_to_10D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_12D_to_11D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_12D_to_11D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_err_z
    
    ! Internal implementation
    pure subroutine norm_12D_to_11D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 12) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_to_11D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_13D_to_12D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_13D_to_12D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_err_z
    
    ! Internal implementation
    pure subroutine norm_13D_to_12D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 13) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_to_12D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_14D_to_13D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_14D_to_13D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_err_z
    
    ! Internal implementation
    pure subroutine norm_14D_to_13D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 14) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_to_13D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_15D_to_14D_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_z(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_z

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_15D_to_14D_err_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_z(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_err_z
    
    ! Internal implementation
    pure subroutine norm_15D_to_14D_z(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 15) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_to_14D_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_2D_to_1D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_2D_to_1D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_2D_to_1D_err_w
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 2) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_2D_to_1D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_3D_to_2D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_3D_to_2D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_3D_to_2D_err_w
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 3) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_3D_to_2D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_4D_to_3D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_4D_to_3D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_4D_to_3D_err_w
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 4) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_4D_to_3D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_5D_to_4D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_5D_to_4D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_5D_to_4D_err_w
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 5) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_5D_to_4D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_6D_to_5D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_6D_to_5D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_6D_to_5D_err_w
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 6) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_6D_to_5D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_7D_to_6D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_7D_to_6D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_7D_to_6D_err_w
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 7) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_7D_to_6D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_8D_to_7D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_8D_to_7D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        call norm_8D_to_7D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_8D_to_7D_err_w
    
    ! Internal implementation
    pure subroutine norm_8D_to_7D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 8) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_8D_to_7D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_9D_to_8D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_9D_to_8D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim))
        
        call norm_9D_to_8D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_9D_to_8D_err_w
    
    ! Internal implementation
    pure subroutine norm_9D_to_8D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 9) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_9D_to_8D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_10D_to_9D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_10D_to_9D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        call norm_10D_to_9D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_10D_to_9D_err_w
    
    ! Internal implementation
    pure subroutine norm_10D_to_9D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 10) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_10D_to_9D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_11D_to_10D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_11D_to_10D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))
        
        call norm_11D_to_10D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_11D_to_10D_err_w
    
    ! Internal implementation
    pure subroutine norm_11D_to_10D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim))  &
            &   
type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 11) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_11D_to_10D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_12D_to_11D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_12D_to_11D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim))
        
        call norm_12D_to_11D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_12D_to_11D_err_w
    
    ! Internal implementation
    pure subroutine norm_12D_to_11D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 12) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_12D_to_11D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_13D_to_12D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_13D_to_12D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        call norm_13D_to_12D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_13D_to_12D_err_w
    
    ! Internal implementation
    pure subroutine norm_13D_to_12D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 13) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_13D_to_12D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_14D_to_13D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_14D_to_13D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim))
        
        call norm_14D_to_13D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_14D_to_13D_err_w
    
    ! Internal implementation
    pure subroutine norm_14D_to_13D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 14) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_14D_to_13D_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function stdlib_linalg_norm_15D_to_14D_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_w(a,order,dim,nrm=nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_w

    ! Function interface with DIM specifier and output error state.
    function stdlib_linalg_norm_15D_to_14D_err_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> Output state return flag.
        type(linalg_state),intent(out) :: err
        !> Norm of the matrix.
        complex(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),size(a,9),&
            & mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),merge(size(a,&
            & 11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),size(a,14),&
            & mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        call norm_15D_to_14D_w(a,order,dim,err,nrm)
            
    end function stdlib_linalg_norm_15D_to_14D_err_w
    
    ! Internal implementation
    pure subroutine norm_15D_to_14D_w(a,order,dim,err,nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer,intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer,intent(in) :: dim
        !> [optional] state return flag. On error if not requested, the code will stop
        type(linalg_state),intent(out),optional :: err
        !> Norm of the matrix.
        complex(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim),merge(size(a,7),size(a,8),mask=7 < dim),merge(size(a,8),&
            & size(a,9),mask=8 < dim),merge(size(a,9),size(a,10),mask=9 < dim),merge(size(a,10),size(a,11),mask=10 < dim),&
            & merge(size(a,11),size(a,12),mask=11 < dim),merge(size(a,12),size(a,13),mask=12 < dim),merge(size(a,13),&
            & size(a,14),mask=13 < dim),merge(size(a,14),size(a,15),mask=14 < dim))
        
        type(linalg_state) :: err_
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if

        if (dim < 1 .or. dim > 15) then
            err_ = linalg_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call linalg_error_handling(err_,err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err_)
        if (err_%error()) then
            call linalg_error_handling(err_,err)
            return
        end if

        select case (order)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(sum(a*conjg(a),dim=dim))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(abs(a)**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine norm_15D_to_14D_w

end module stdlib_linalg_norms
