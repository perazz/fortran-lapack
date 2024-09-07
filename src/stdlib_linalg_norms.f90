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
        !> Scalar norms
        module procedure stdlib_linalg_norm_1D_s
        module procedure stdlib_linalg_norm_2D_s
        module procedure stdlib_linalg_norm_3D_s
        module procedure stdlib_linalg_norm_4D_s
        module procedure stdlib_linalg_norm_5D_s
        module procedure stdlib_linalg_norm_6D_s
        module procedure stdlib_linalg_norm_7D_s
        module procedure stdlib_linalg_norm_2D_to_1D_s
        module procedure stdlib_linalg_norm_3D_to_2D_s
        module procedure stdlib_linalg_norm_4D_to_3D_s
        module procedure stdlib_linalg_norm_5D_to_4D_s
        module procedure stdlib_linalg_norm_6D_to_5D_s
        module procedure stdlib_linalg_norm_7D_to_6D_s
        module procedure stdlib_linalg_norm_1D_d
        module procedure stdlib_linalg_norm_2D_d
        module procedure stdlib_linalg_norm_3D_d
        module procedure stdlib_linalg_norm_4D_d
        module procedure stdlib_linalg_norm_5D_d
        module procedure stdlib_linalg_norm_6D_d
        module procedure stdlib_linalg_norm_7D_d
        module procedure stdlib_linalg_norm_2D_to_1D_d
        module procedure stdlib_linalg_norm_3D_to_2D_d
        module procedure stdlib_linalg_norm_4D_to_3D_d
        module procedure stdlib_linalg_norm_5D_to_4D_d
        module procedure stdlib_linalg_norm_6D_to_5D_d
        module procedure stdlib_linalg_norm_7D_to_6D_d
        module procedure stdlib_linalg_norm_1D_q
        module procedure stdlib_linalg_norm_2D_q
        module procedure stdlib_linalg_norm_3D_q
        module procedure stdlib_linalg_norm_4D_q
        module procedure stdlib_linalg_norm_5D_q
        module procedure stdlib_linalg_norm_6D_q
        module procedure stdlib_linalg_norm_7D_q
        module procedure stdlib_linalg_norm_2D_to_1D_q
        module procedure stdlib_linalg_norm_3D_to_2D_q
        module procedure stdlib_linalg_norm_4D_to_3D_q
        module procedure stdlib_linalg_norm_5D_to_4D_q
        module procedure stdlib_linalg_norm_6D_to_5D_q
        module procedure stdlib_linalg_norm_7D_to_6D_q
        module procedure stdlib_linalg_norm_1D_c
        module procedure stdlib_linalg_norm_2D_c
        module procedure stdlib_linalg_norm_3D_c
        module procedure stdlib_linalg_norm_4D_c
        module procedure stdlib_linalg_norm_5D_c
        module procedure stdlib_linalg_norm_6D_c
        module procedure stdlib_linalg_norm_7D_c
        module procedure stdlib_linalg_norm_2D_to_1D_c
        module procedure stdlib_linalg_norm_3D_to_2D_c
        module procedure stdlib_linalg_norm_4D_to_3D_c
        module procedure stdlib_linalg_norm_5D_to_4D_c
        module procedure stdlib_linalg_norm_6D_to_5D_c
        module procedure stdlib_linalg_norm_7D_to_6D_c
        module procedure stdlib_linalg_norm_1D_z
        module procedure stdlib_linalg_norm_2D_z
        module procedure stdlib_linalg_norm_3D_z
        module procedure stdlib_linalg_norm_4D_z
        module procedure stdlib_linalg_norm_5D_z
        module procedure stdlib_linalg_norm_6D_z
        module procedure stdlib_linalg_norm_7D_z
        module procedure stdlib_linalg_norm_2D_to_1D_z
        module procedure stdlib_linalg_norm_3D_to_2D_z
        module procedure stdlib_linalg_norm_4D_to_3D_z
        module procedure stdlib_linalg_norm_5D_to_4D_z
        module procedure stdlib_linalg_norm_6D_to_5D_z
        module procedure stdlib_linalg_norm_7D_to_6D_z
        module procedure stdlib_linalg_norm_1D_w
        module procedure stdlib_linalg_norm_2D_w
        module procedure stdlib_linalg_norm_3D_w
        module procedure stdlib_linalg_norm_4D_w
        module procedure stdlib_linalg_norm_5D_w
        module procedure stdlib_linalg_norm_6D_w
        module procedure stdlib_linalg_norm_7D_w
        module procedure stdlib_linalg_norm_2D_to_1D_w
        module procedure stdlib_linalg_norm_3D_to_2D_w
        module procedure stdlib_linalg_norm_4D_to_3D_w
        module procedure stdlib_linalg_norm_5D_to_4D_w
        module procedure stdlib_linalg_norm_6D_to_5D_w
        module procedure stdlib_linalg_norm_7D_to_6D_w
     end interface norm
     
     interface parse_norm_type
        module procedure parse_norm_type_integer
        !module procedure parse_norm_type_character
     end interface

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
     
    !==============================================
    ! Norms : any rank to 0D
    !==============================================

    pure subroutine stdlib_linalg_norm_1D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_1D_s

    pure subroutine stdlib_linalg_norm_2D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_s

    pure subroutine stdlib_linalg_norm_3D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_s

    pure subroutine stdlib_linalg_norm_4D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_s

    pure subroutine stdlib_linalg_norm_5D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_s

    pure subroutine stdlib_linalg_norm_6D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_s

    pure subroutine stdlib_linalg_norm_7D_s(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_s

    pure subroutine stdlib_linalg_norm_1D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_1D_d

    pure subroutine stdlib_linalg_norm_2D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_d

    pure subroutine stdlib_linalg_norm_3D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_d

    pure subroutine stdlib_linalg_norm_4D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_d

    pure subroutine stdlib_linalg_norm_5D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_d

    pure subroutine stdlib_linalg_norm_6D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_d

    pure subroutine stdlib_linalg_norm_7D_d(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_d

    pure subroutine stdlib_linalg_norm_1D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_1D_q

    pure subroutine stdlib_linalg_norm_2D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_q

    pure subroutine stdlib_linalg_norm_3D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_q

    pure subroutine stdlib_linalg_norm_4D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_q

    pure subroutine stdlib_linalg_norm_5D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_q

    pure subroutine stdlib_linalg_norm_6D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_q

    pure subroutine stdlib_linalg_norm_7D_q(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_q

    pure subroutine stdlib_linalg_norm_1D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_1D_c

    pure subroutine stdlib_linalg_norm_2D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_c

    pure subroutine stdlib_linalg_norm_3D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_c

    pure subroutine stdlib_linalg_norm_4D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_c

    pure subroutine stdlib_linalg_norm_5D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_c

    pure subroutine stdlib_linalg_norm_6D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_c

    pure subroutine stdlib_linalg_norm_7D_c(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_c

    pure subroutine stdlib_linalg_norm_1D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_1D_z

    pure subroutine stdlib_linalg_norm_2D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_z

    pure subroutine stdlib_linalg_norm_3D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_z

    pure subroutine stdlib_linalg_norm_4D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_z

    pure subroutine stdlib_linalg_norm_5D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_z

    pure subroutine stdlib_linalg_norm_6D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_z

    pure subroutine stdlib_linalg_norm_7D_z(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_z

    pure subroutine stdlib_linalg_norm_1D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_1D_w

    pure subroutine stdlib_linalg_norm_2D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_w

    pure subroutine stdlib_linalg_norm_3D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_w

    pure subroutine stdlib_linalg_norm_4D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_w

    pure subroutine stdlib_linalg_norm_5D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_w

    pure subroutine stdlib_linalg_norm_6D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_w

    pure subroutine stdlib_linalg_norm_7D_w(a,order,err,nrm)
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
            case (NORM_TWO)
                nrm = nrm2(sze,a,1_ilp)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_w

    !==============================================
    ! Norms : any rank to rank-1, for rank > 2
    !==============================================

    pure subroutine stdlib_linalg_norm_2D_to_1D_s(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_to_1D_s

    pure subroutine stdlib_linalg_norm_3D_to_2D_s(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_to_2D_s

    pure subroutine stdlib_linalg_norm_4D_to_3D_s(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_to_3D_s

    pure subroutine stdlib_linalg_norm_5D_to_4D_s(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_to_4D_s

    pure subroutine stdlib_linalg_norm_6D_to_5D_s(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_to_5D_s

    pure subroutine stdlib_linalg_norm_7D_to_6D_s(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_to_6D_s

    pure subroutine stdlib_linalg_norm_2D_to_1D_d(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_to_1D_d

    pure subroutine stdlib_linalg_norm_3D_to_2D_d(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_to_2D_d

    pure subroutine stdlib_linalg_norm_4D_to_3D_d(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_to_3D_d

    pure subroutine stdlib_linalg_norm_5D_to_4D_d(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_to_4D_d

    pure subroutine stdlib_linalg_norm_6D_to_5D_d(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_to_5D_d

    pure subroutine stdlib_linalg_norm_7D_to_6D_d(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_to_6D_d

    pure subroutine stdlib_linalg_norm_2D_to_1D_q(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_to_1D_q

    pure subroutine stdlib_linalg_norm_3D_to_2D_q(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_to_2D_q

    pure subroutine stdlib_linalg_norm_4D_to_3D_q(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_to_3D_q

    pure subroutine stdlib_linalg_norm_5D_to_4D_q(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_to_4D_q

    pure subroutine stdlib_linalg_norm_6D_to_5D_q(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_to_5D_q

    pure subroutine stdlib_linalg_norm_7D_to_6D_q(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_to_6D_q

    pure subroutine stdlib_linalg_norm_2D_to_1D_c(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_to_1D_c

    pure subroutine stdlib_linalg_norm_3D_to_2D_c(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_to_2D_c

    pure subroutine stdlib_linalg_norm_4D_to_3D_c(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_to_3D_c

    pure subroutine stdlib_linalg_norm_5D_to_4D_c(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_to_4D_c

    pure subroutine stdlib_linalg_norm_6D_to_5D_c(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_to_5D_c

    pure subroutine stdlib_linalg_norm_7D_to_6D_c(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_sp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_to_6D_c

    pure subroutine stdlib_linalg_norm_2D_to_1D_z(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_to_1D_z

    pure subroutine stdlib_linalg_norm_3D_to_2D_z(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_to_2D_z

    pure subroutine stdlib_linalg_norm_4D_to_3D_z(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_to_3D_z

    pure subroutine stdlib_linalg_norm_5D_to_4D_z(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_to_4D_z

    pure subroutine stdlib_linalg_norm_6D_to_5D_z(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_to_5D_z

    pure subroutine stdlib_linalg_norm_7D_to_6D_z(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_dp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_to_6D_z

    pure subroutine stdlib_linalg_norm_2D_to_1D_w(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_2D_to_1D_w

    pure subroutine stdlib_linalg_norm_3D_to_2D_w(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_3D_to_2D_w

    pure subroutine stdlib_linalg_norm_4D_to_3D_w(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_4D_to_3D_w

    pure subroutine stdlib_linalg_norm_5D_to_4D_w(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_5D_to_4D_w

    pure subroutine stdlib_linalg_norm_6D_to_5D_w(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_6D_to_5D_w

    pure subroutine stdlib_linalg_norm_7D_to_6D_w(a,order,dim,err,nrm)
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
            case (NORM_TWO)
                nrm = sqrt(sum(a**2,dim=dim))
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_P)
                rorder = 1.0_qp/order
                nrm = sum(a**order,dim=dim)**rorder
            case default
                err_ = linalg_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call linalg_error_handling(err_,err)
        end select
        
    end subroutine stdlib_linalg_norm_7D_to_6D_w

end module stdlib_linalg_norms
