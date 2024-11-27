! Test inverse matrix
module test_linalg_pseudoinverse
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> Matrix inversion tests
    subroutine test_pseudoinverse_matrix(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_s_eye_pseudoinverse(error)
        if (error) return
        call test_d_eye_pseudoinverse(error)
        if (error) return
        call test_q_eye_pseudoinverse(error)
        if (error) return
        call test_s_random_pseudoinverse(error)
        if (error) return
        call test_d_random_pseudoinverse(error)
        if (error) return
        call test_q_random_pseudoinverse(error)
        if (error) return
        call test_c_random_pseudoinverse(error)
        if (error) return
        call test_z_random_pseudoinverse(error)
        if (error) return
        call test_w_random_pseudoinverse(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Pseudo-Inverse matrix tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_pseudoinverse_matrix

    !> Invert identity matrix
    subroutine test_s_eye_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))

        real(sp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_sp,0.0_sp,i == j)
        end do

        !> Invert function
        inva = pinv(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tol)
        if (error) return

        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tol)
        
        !> Operator
        inva = .pinv.a
        error = .not. all(abs(a - inva) < tiny(0.0_sp))

    end subroutine test_s_eye_pseudoinverse

    subroutine test_d_eye_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(dp),parameter :: tol = sqrt(epsilon(0.0_sp))

        real(dp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_dp,0.0_dp,i == j)
        end do

        !> Invert function
        inva = pinv(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tol)
        if (error) return

        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tol)
        
        !> Operator
        inva = .pinv.a
        error = .not. all(abs(a - inva) < tiny(0.0_dp))

    end subroutine test_d_eye_pseudoinverse

    subroutine test_q_eye_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))

        real(qp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_qp,0.0_qp,i == j)
        end do

        !> Invert function
        inva = pinv(a,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tol)
        if (error) return

        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        error = state%error() .or. .not. all(abs(a - inva) < tol)
        
        !> Operator
        inva = .pinv.a
        error = .not. all(abs(a - inva) < tiny(0.0_qp))

    end subroutine test_q_eye_pseudoinverse

    !> Invert random matrix
    subroutine test_s_random_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 7,n = 15
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        error = state%error(); if (error) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        print *, 'failed=',failed
        error = failed > 0; if (error) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        print *, 'failed=',failed
        error = failed > 0; if (error) return

    end subroutine test_s_random_pseudoinverse

    subroutine test_d_random_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 7,n = 15
        real(dp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(dp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        error = state%error(); if (error) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        error = failed > 0; if (error) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        error = failed > 0; if (error) return

    end subroutine test_d_random_pseudoinverse

    subroutine test_q_random_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 7,n = 15
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        error = state%error(); if (error) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        error = failed > 0; if (error) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        error = failed > 0; if (error) return

    end subroutine test_q_random_pseudoinverse

    subroutine test_c_random_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 7,n = 15
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),inva(n,m)
        real(sp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=sp)
        
        inva = pinv(a,err=state)
        error = state%error(); if (error) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        error = failed > 0; if (error) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        error = failed > 0; if (error) return

    end subroutine test_c_random_pseudoinverse

    subroutine test_z_random_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 7,n = 15
        real(dp),parameter :: tol = sqrt(epsilon(0.0_sp))
        complex(dp) :: a(m,n),inva(n,m)
        real(dp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=dp)
        
        inva = pinv(a,err=state)
        error = state%error(); if (error) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        error = failed > 0; if (error) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        error = failed > 0; if (error) return

    end subroutine test_z_random_pseudoinverse

    subroutine test_w_random_pseudoinverse(error)
        logical,intent(out) :: error

        type(linalg_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 7,n = 15
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),inva(n,m)
        real(qp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=qp)
        
        inva = pinv(a,err=state)
        error = state%error(); if (error) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        error = failed > 0; if (error) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        error = failed > 0; if (error) return

    end subroutine test_w_random_pseudoinverse

end module test_linalg_pseudoinverse

