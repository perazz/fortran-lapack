! Test inverse matrix
module test_linalg_inverse
    use la_interface

    implicit none(type,external)

    contains

    !> Matrix inversion tests
    subroutine test_inverse_matrix(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

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

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Inverse matrix tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_inverse_matrix

    !> Invert identity matrix
    subroutine test_s_eye_inverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 250_ilp

        real(sp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_sp,0.0_sp,i == j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tiny(0.0_sp))
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tiny(0.0_sp))

    end subroutine test_s_eye_inverse

    subroutine test_d_eye_inverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 250_ilp

        real(dp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_dp,0.0_dp,i == j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tiny(0.0_dp))
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tiny(0.0_dp))

    end subroutine test_d_eye_inverse

    subroutine test_q_eye_inverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 250_ilp

        real(qp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_qp,0.0_qp,i == j)
        end do

        !> Invert function
        inva = inv(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tiny(0.0_qp))
        if (error) return

        !> Inverse subroutine
        call invert(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tiny(0.0_qp))

    end subroutine test_q_eye_inverse

    !> Invert identity matrix
    subroutine test_c_eye_inverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 250_ilp

        complex(sp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_sp,1.0_sp), (0.0_sp,0.0_sp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        error = state%error() .or. .not. failed == 0
        if (error) return

        !> Inverse subroutine
        call invert(copya,err=state)

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        error = state%error() .or. .not. failed == 0

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(sp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < tiny(0.0_sp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_sp,0.0_sp)) < tiny(0.0_sp)
               end if
           end function is_diagonal_inverse

    end subroutine test_c_eye_inverse

    subroutine test_z_eye_inverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 250_ilp

        complex(dp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_dp,1.0_dp), (0.0_dp,0.0_dp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        error = state%error() .or. .not. failed == 0
        if (error) return

        !> Inverse subroutine
        call invert(copya,err=state)

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        error = state%error() .or. .not. failed == 0

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(dp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < tiny(0.0_dp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_dp,0.0_dp)) < tiny(0.0_dp)
               end if
           end function is_diagonal_inverse

    end subroutine test_z_eye_inverse

    subroutine test_w_eye_inverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 250_ilp

        complex(qp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_qp,1.0_qp), (0.0_qp,0.0_qp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        error = state%error() .or. .not. failed == 0
        if (error) return

        !> Inverse subroutine
        call invert(copya,err=state)

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        error = state%error() .or. .not. failed == 0

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(qp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < tiny(0.0_qp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_qp,0.0_qp)) < tiny(0.0_qp)
               end if
           end function is_diagonal_inverse

    end subroutine test_w_eye_inverse

end module test_linalg_inverse

