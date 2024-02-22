
module test_linalg_eye
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> Identity matrix tests
    subroutine test_eye(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_s_eye_allocation(error)
        if (error) return
        call test_d_eye_allocation(error)
        if (error) return
        call test_q_eye_allocation(error)
        if (error) return
        call test_c_eye_allocation(error)
        if (error) return
        call test_z_eye_allocation(error)
        if (error) return
        call test_w_eye_allocation(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Identity matrix tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_eye

    !> Determinant of identity matrix
    subroutine test_s_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(sp),allocatable :: a(:,:)
        real(sp) :: dummy

        !> Should be error
        a = eye(-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,dtype=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

    end subroutine test_s_eye_allocation

    subroutine test_d_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(dp),allocatable :: a(:,:)
        real(dp) :: dummy

        !> Should be error
        a = eye(-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,dtype=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

    end subroutine test_d_eye_allocation

    subroutine test_q_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(qp),allocatable :: a(:,:)
        real(qp) :: dummy

        !> Should be error
        a = eye(-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,dtype=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

    end subroutine test_q_eye_allocation

    subroutine test_c_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(sp),allocatable :: a(:,:)
        complex(sp) :: dummy

        !> Should be error
        a = eye(-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,dtype=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

    end subroutine test_c_eye_allocation

    subroutine test_z_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(dp),allocatable :: a(:,:)
        complex(dp) :: dummy

        !> Should be error
        a = eye(-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,dtype=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

    end subroutine test_z_eye_allocation

    subroutine test_w_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(qp),allocatable :: a(:,:)
        complex(qp) :: dummy

        !> Should be error
        a = eye(-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,dtype=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,dtype=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,dtype=dummy,err=state)
        error = state%error() .or. nint(sum(a),kind=ilp) /= 5
        if (error) return

    end subroutine test_w_eye_allocation

end module test_linalg_eye

