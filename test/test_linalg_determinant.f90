! Test matrix determinant
module test_linalg_determinant
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> Matrix inversion tests
    subroutine test_matrix_determinant(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_s_eye_determinant(error)
        if (error) return
        call test_s_eye_multiple(error)
        if (error) return
        call test_d_eye_determinant(error)
        if (error) return
        call test_d_eye_multiple(error)
        if (error) return
        call test_q_eye_determinant(error)
        if (error) return
        call test_q_eye_multiple(error)
        if (error) return
!        call test_c_eye_inverse(error)
!        if (error) return
!        call test_z_eye_inverse(error)
!        if (error) return
!        call test_w_eye_inverse(error)
!        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Determinant tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_matrix_determinant

    !> Determinant of identity matrix
    subroutine test_s_eye_determinant(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 128_ilp

        real(sp) :: a(n,n),deta

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_sp,0.0_sp,i == j)
        end do

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_sp) < tiny(0.0_sp)
        if (error) return

    end subroutine test_s_eye_determinant

    subroutine test_s_eye_multiple(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(sp),parameter :: coef = 0.0001_sp
        integer(ilp) :: i,j
        real(sp) :: a(n,n),deta

        !> Multiply eye by a very small number
        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(coef,0.0_sp,i == j)
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < tiny(0.0_sp)
        if (error) return

    end subroutine test_s_eye_multiple

    subroutine test_d_eye_determinant(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 128_ilp

        real(dp) :: a(n,n),deta

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_dp,0.0_dp,i == j)
        end do

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_dp) < tiny(0.0_dp)
        if (error) return

    end subroutine test_d_eye_determinant

    subroutine test_d_eye_multiple(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(dp),parameter :: coef = 0.0001_dp
        integer(ilp) :: i,j
        real(dp) :: a(n,n),deta

        !> Multiply eye by a very small number
        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(coef,0.0_dp,i == j)
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < tiny(0.0_dp)
        if (error) return

    end subroutine test_d_eye_multiple

    subroutine test_q_eye_determinant(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 128_ilp

        real(qp) :: a(n,n),deta

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_qp,0.0_qp,i == j)
        end do

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_qp) < tiny(0.0_qp)
        if (error) return

    end subroutine test_q_eye_determinant

    subroutine test_q_eye_multiple(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(qp),parameter :: coef = 0.0001_qp
        integer(ilp) :: i,j
        real(qp) :: a(n,n),deta

        !> Multiply eye by a very small number
        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(coef,0.0_qp,i == j)
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < tiny(0.0_qp)
        if (error) return

    end subroutine test_q_eye_multiple

    !> Determinant of identity matrix
    subroutine test_c_eye_determinant(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 250_ilp

        complex(sp) :: a(n,n),copya(n,n),inva(n,n)

        ! Get eye matrix
        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_sp,1.0_sp), (0.0_sp,0.0_sp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = det(a,err=state)

        error = state%error() .or. .not. failed == 0
        if (error) return

        !> Inverse subroutine
        call invert(copya,err=state)
!
!        failed = 0
!        do i=1,n
!            do j=1,n
!                if (.not.is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed+1
!            end do
!        end do

        error = state%error() .or. .not. failed == 0

    end subroutine test_c_eye_determinant

    subroutine test_z_eye_determinant(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 250_ilp

        complex(dp) :: a(n,n),copya(n,n),inva(n,n)

        ! Get eye matrix
        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_dp,1.0_dp), (0.0_dp,0.0_dp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = det(a,err=state)

        error = state%error() .or. .not. failed == 0
        if (error) return

        !> Inverse subroutine
        call invert(copya,err=state)
!
!        failed = 0
!        do i=1,n
!            do j=1,n
!                if (.not.is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed+1
!            end do
!        end do

        error = state%error() .or. .not. failed == 0

    end subroutine test_z_eye_determinant

    subroutine test_w_eye_determinant(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 250_ilp

        complex(qp) :: a(n,n),copya(n,n),inva(n,n)

        ! Get eye matrix
        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_qp,1.0_qp), (0.0_qp,0.0_qp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = det(a,err=state)

        error = state%error() .or. .not. failed == 0
        if (error) return

        !> Inverse subroutine
        call invert(copya,err=state)
!
!        failed = 0
!        do i=1,n
!            do j=1,n
!                if (.not.is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed+1
!            end do
!        end do

        error = state%error() .or. .not. failed == 0

    end subroutine test_w_eye_determinant

end module test_linalg_determinant

