! Test inverse matrix
module test_linalg_inverse
    use stdlib_linalg_interface

    implicit none (type,external)

    contains

    !> Matrix inversion tests
    subroutine test_inverse_matrix(error)
        logical, intent(out) :: error

        call test_s_eye_inverse(error)
        if (error) return
        call test_d_eye_inverse(error)
        if (error) return
        call test_q_eye_inverse(error)
        if (error) return
        call test_c_eye_inverse(error)
        if (error) return
        call test_z_eye_inverse(error)
        if (error) return
        call test_w_eye_inverse(error)
        if (error) return

    end subroutine test_inverse_matrix

    !> Invert identity matrix
    subroutine test_s_eye_inverse(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp), parameter :: n = 500_ilp

        real(sp) :: a(n,n),inva(n,n)

        do concurrent (i=1:n,j=1:n)
          a(i,j) = merge(1.0_sp,0.0_sp,i==j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not.all(abs(a-inva)==0.0_sp)
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not.all(a==inva)

    end subroutine test_s_eye_inverse

    subroutine test_d_eye_inverse(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp), parameter :: n = 500_ilp

        real(dp) :: a(n,n),inva(n,n)

        do concurrent (i=1:n,j=1:n)
          a(i,j) = merge(1.0_dp,0.0_dp,i==j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not.all(abs(a-inva)==0.0_dp)
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not.all(a==inva)

    end subroutine test_d_eye_inverse

    subroutine test_q_eye_inverse(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp), parameter :: n = 500_ilp

        real(qp) :: a(n,n),inva(n,n)

        do concurrent (i=1:n,j=1:n)
          a(i,j) = merge(1.0_qp,0.0_qp,i==j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not.all(abs(a-inva)==0.0_qp)
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not.all(a==inva)

    end subroutine test_q_eye_inverse


    !> Invert identity matrix
    subroutine test_c_eye_inverse(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp), parameter :: n = 500_ilp

        complex(sp) :: a(n,n),inva(n,n)

        do concurrent (i=1:n,j=1:n)
          a(i,j) = merge((1.0_sp,1.0_sp),(0.0_sp,0.0_sp),i==j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not.all(abs(a-inva)==0.0_sp)
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not.all(a==inva)

    end subroutine test_c_eye_inverse

    subroutine test_z_eye_inverse(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp), parameter :: n = 500_ilp

        complex(dp) :: a(n,n),inva(n,n)

        do concurrent (i=1:n,j=1:n)
          a(i,j) = merge((1.0_dp,1.0_dp),(0.0_dp,0.0_dp),i==j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not.all(abs(a-inva)==0.0_dp)
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not.all(a==inva)

    end subroutine test_z_eye_inverse

    subroutine test_w_eye_inverse(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp), parameter :: n = 500_ilp

        complex(qp) :: a(n,n),inva(n,n)

        do concurrent (i=1:n,j=1:n)
          a(i,j) = merge((1.0_qp,1.0_qp),(0.0_qp,0.0_qp),i==j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not.all(abs(a-inva)==0.0_qp)
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not.all(a==inva)

    end subroutine test_w_eye_inverse


end module test_linalg_inverse


