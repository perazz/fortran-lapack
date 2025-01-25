module test_linalg_eye
    use la_interface

    implicit none(type,external)

    contains

    !> Identity matrix tests
    subroutine test_eye(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_s_eye_allocation(error)
        call test_s_diag_scalar(error)
        call test_s_diag_array(error)
        if (error) return
        call test_d_eye_allocation(error)
        call test_d_diag_scalar(error)
        call test_d_diag_array(error)
        if (error) return
        call test_q_eye_allocation(error)
        call test_q_diag_scalar(error)
        call test_q_diag_array(error)
        if (error) return
        call test_c_eye_allocation(error)
        call test_c_diag_scalar(error)
        call test_c_diag_array(error)
        if (error) return
        call test_z_eye_allocation(error)
        call test_z_diag_scalar(error)
        call test_z_diag_array(error)
        if (error) return
        call test_w_eye_allocation(error)
        call test_w_diag_scalar(error)
        call test_w_diag_array(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Identity matrix tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_eye

    !> Identity matrix: test allocation
    subroutine test_s_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(sp),allocatable :: a(:,:)
        real(sp) :: dummy

        !> Should be error
        a = eye(-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,mold=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 5
        if (error) return

    end subroutine test_s_eye_allocation

    !> Diagonal matrix from scalar
    subroutine test_s_diag_scalar(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(sp),allocatable :: a(:,:)
        real(sp) :: dummy

        dummy = 2.0_sp

        !> Should be error
        a = diag(-1,source=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = diag(0,source=dummy,err=state)
        error = state%error() .or. any(shape(a) /= 0)
        if (error) return

        !> Test identity values
        a = diag(5,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 10

        a = diag(12,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 24
        if (error) return

    end subroutine test_s_diag_scalar

    !> Diagonal matrix from array
    subroutine test_s_diag_array(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        integer(ilp),parameter :: array_sizes(*) = [1,2,5,10,20,50,100]
        real(sp),allocatable :: a(:,:),darr(:)

        do i = 1,size(array_sizes)

            allocate (darr(array_sizes(i)))
            darr = 2.0_sp

            !> Test with several sizes
            a = diag(source=darr,err=state)
            error = state%error() .or. any(shape(a) /= size(darr)) .or. sum(a) /= sum(darr)
            if (error) return

            deallocate (darr)

        end do

    end subroutine test_s_diag_array

    !> Identity matrix: test allocation
    subroutine test_d_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(dp),allocatable :: a(:,:)
        real(dp) :: dummy

        !> Should be error
        a = eye(-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,mold=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 5
        if (error) return

    end subroutine test_d_eye_allocation

    !> Diagonal matrix from scalar
    subroutine test_d_diag_scalar(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(dp),allocatable :: a(:,:)
        real(dp) :: dummy

        dummy = 2.0_dp

        !> Should be error
        a = diag(-1,source=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = diag(0,source=dummy,err=state)
        error = state%error() .or. any(shape(a) /= 0)
        if (error) return

        !> Test identity values
        a = diag(5,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 10

        a = diag(12,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 24
        if (error) return

    end subroutine test_d_diag_scalar

    !> Diagonal matrix from array
    subroutine test_d_diag_array(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        integer(ilp),parameter :: array_sizes(*) = [1,2,5,10,20,50,100]
        real(dp),allocatable :: a(:,:),darr(:)

        do i = 1,size(array_sizes)

            allocate (darr(array_sizes(i)))
            darr = 2.0_dp

            !> Test with several sizes
            a = diag(source=darr,err=state)
            error = state%error() .or. any(shape(a) /= size(darr)) .or. sum(a) /= sum(darr)
            if (error) return

            deallocate (darr)

        end do

    end subroutine test_d_diag_array

    !> Identity matrix: test allocation
    subroutine test_q_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(qp),allocatable :: a(:,:)
        real(qp) :: dummy

        !> Should be error
        a = eye(-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,mold=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 5
        if (error) return

    end subroutine test_q_eye_allocation

    !> Diagonal matrix from scalar
    subroutine test_q_diag_scalar(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        real(qp),allocatable :: a(:,:)
        real(qp) :: dummy

        dummy = 2.0_qp

        !> Should be error
        a = diag(-1,source=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = diag(0,source=dummy,err=state)
        error = state%error() .or. any(shape(a) /= 0)
        if (error) return

        !> Test identity values
        a = diag(5,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 10

        a = diag(12,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 24
        if (error) return

    end subroutine test_q_diag_scalar

    !> Diagonal matrix from array
    subroutine test_q_diag_array(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        integer(ilp),parameter :: array_sizes(*) = [1,2,5,10,20,50,100]
        real(qp),allocatable :: a(:,:),darr(:)

        do i = 1,size(array_sizes)

            allocate (darr(array_sizes(i)))
            darr = 2.0_qp

            !> Test with several sizes
            a = diag(source=darr,err=state)
            error = state%error() .or. any(shape(a) /= size(darr)) .or. sum(a) /= sum(darr)
            if (error) return

            deallocate (darr)

        end do

    end subroutine test_q_diag_array

    !> Identity matrix: test allocation
    subroutine test_c_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(sp),allocatable :: a(:,:)
        complex(sp) :: dummy

        !> Should be error
        a = eye(-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,mold=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 5
        if (error) return

    end subroutine test_c_eye_allocation

    !> Diagonal matrix from scalar
    subroutine test_c_diag_scalar(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(sp),allocatable :: a(:,:)
        complex(sp) :: dummy

        dummy = 2.0_sp

        !> Should be error
        a = diag(-1,source=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = diag(0,source=dummy,err=state)
        error = state%error() .or. any(shape(a) /= 0)
        if (error) return

        !> Test identity values
        a = diag(5,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 10

        a = diag(12,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=sp),kind=ilp) /= 24
        if (error) return

    end subroutine test_c_diag_scalar

    !> Diagonal matrix from array
    subroutine test_c_diag_array(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        integer(ilp),parameter :: array_sizes(*) = [1,2,5,10,20,50,100]
        complex(sp),allocatable :: a(:,:),darr(:)

        do i = 1,size(array_sizes)

            allocate (darr(array_sizes(i)))
            darr = 2.0_sp

            !> Test with several sizes
            a = diag(source=darr,err=state)
            error = state%error() .or. any(shape(a) /= size(darr)) .or. sum(a) /= sum(darr)
            if (error) return

            deallocate (darr)

        end do

    end subroutine test_c_diag_array

    !> Identity matrix: test allocation
    subroutine test_z_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(dp),allocatable :: a(:,:)
        complex(dp) :: dummy

        !> Should be error
        a = eye(-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,mold=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 5
        if (error) return

    end subroutine test_z_eye_allocation

    !> Diagonal matrix from scalar
    subroutine test_z_diag_scalar(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(dp),allocatable :: a(:,:)
        complex(dp) :: dummy

        dummy = 2.0_dp

        !> Should be error
        a = diag(-1,source=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = diag(0,source=dummy,err=state)
        error = state%error() .or. any(shape(a) /= 0)
        if (error) return

        !> Test identity values
        a = diag(5,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 10

        a = diag(12,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=dp),kind=ilp) /= 24
        if (error) return

    end subroutine test_z_diag_scalar

    !> Diagonal matrix from array
    subroutine test_z_diag_array(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        integer(ilp),parameter :: array_sizes(*) = [1,2,5,10,20,50,100]
        complex(dp),allocatable :: a(:,:),darr(:)

        do i = 1,size(array_sizes)

            allocate (darr(array_sizes(i)))
            darr = 2.0_dp

            !> Test with several sizes
            a = diag(source=darr,err=state)
            error = state%error() .or. any(shape(a) /= size(darr)) .or. sum(a) /= sum(darr)
            if (error) return

            deallocate (darr)

        end do

    end subroutine test_z_diag_array

    !> Identity matrix: test allocation
    subroutine test_w_eye_allocation(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(qp),allocatable :: a(:,:)
        complex(qp) :: dummy

        !> Should be error
        a = eye(-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be error
        a = eye(5,-1,mold=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = eye(0,5,mold=dummy,err=state)
        error = state%error() .or. any(shape(a) /= [0,5])
        if (error) return

        !> Test identity values
        a = eye(5,10,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 5
        if (error) return

        a = eye(10,5,mold=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 5
        if (error) return

    end subroutine test_w_eye_allocation

    !> Diagonal matrix from scalar
    subroutine test_w_diag_scalar(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        complex(qp),allocatable :: a(:,:)
        complex(qp) :: dummy

        dummy = 2.0_qp

        !> Should be error
        a = diag(-1,source=dummy,err=state)
        error = .not. state%error()
        if (error) return

        !> Should be ok
        a = diag(0,source=dummy,err=state)
        error = state%error() .or. any(shape(a) /= 0)
        if (error) return

        !> Test identity values
        a = diag(5,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 10

        a = diag(12,source=dummy,err=state)
        error = state%error() .or. nint(real(sum(a),kind=qp),kind=ilp) /= 24
        if (error) return

    end subroutine test_w_diag_scalar

    !> Diagonal matrix from array
    subroutine test_w_diag_array(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i

        integer(ilp),parameter :: array_sizes(*) = [1,2,5,10,20,50,100]
        complex(qp),allocatable :: a(:,:),darr(:)

        do i = 1,size(array_sizes)

            allocate (darr(array_sizes(i)))
            darr = 2.0_qp

            !> Test with several sizes
            a = diag(source=darr,err=state)
            error = state%error() .or. any(shape(a) /= size(darr)) .or. sum(a) /= sum(darr)
            if (error) return

            deallocate (darr)

        end do

    end subroutine test_w_diag_array

end module test_linalg_eye

