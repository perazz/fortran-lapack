! Test Moore-Penrose pseudo matrix inverse
module test_linalg_pseudoinverse
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use linear_algebra
    use la_constants

    implicit none(type,external)
    private
    
    public :: test_pseudoinverse_matrix

    contains

    !> Matrix pseudo-inversion tests
    subroutine test_pseudoinverse_matrix(tests)
        !> Collertion of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)
        
        allocate (tests(0))

        call add_test(tests,new_unittest("s_eye_pseudoinverse",test_s_eye_pseudoinverse))
        call add_test(tests,new_unittest("d_eye_pseudoinverse",test_d_eye_pseudoinverse))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("x_eye_pseudoinverse",test_x_eye_pseudoinverse))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("q_eye_pseudoinverse",test_q_eye_pseudoinverse))
#endif
        call add_test(tests,new_unittest("s_square_pseudoinverse",test_s_square_pseudoinverse))
        call add_test(tests,new_unittest("s_tall_pseudoinverse",test_s_tall_pseudoinverse))
        call add_test(tests,new_unittest("s_wide_pseudoinverse",test_s_wide_pseudoinverse))
        call add_test(tests,new_unittest("s_singular_pseudoinverse",test_s_singular_pseudoinverse))
        call add_test(tests,new_unittest("d_square_pseudoinverse",test_d_square_pseudoinverse))
        call add_test(tests,new_unittest("d_tall_pseudoinverse",test_d_tall_pseudoinverse))
        call add_test(tests,new_unittest("d_wide_pseudoinverse",test_d_wide_pseudoinverse))
        call add_test(tests,new_unittest("d_singular_pseudoinverse",test_d_singular_pseudoinverse))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("x_square_pseudoinverse",test_x_square_pseudoinverse))
        call add_test(tests,new_unittest("x_tall_pseudoinverse",test_x_tall_pseudoinverse))
        call add_test(tests,new_unittest("x_wide_pseudoinverse",test_x_wide_pseudoinverse))
        call add_test(tests,new_unittest("x_singular_pseudoinverse",test_x_singular_pseudoinverse))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("q_square_pseudoinverse",test_q_square_pseudoinverse))
        call add_test(tests,new_unittest("q_tall_pseudoinverse",test_q_tall_pseudoinverse))
        call add_test(tests,new_unittest("q_wide_pseudoinverse",test_q_wide_pseudoinverse))
        call add_test(tests,new_unittest("q_singular_pseudoinverse",test_q_singular_pseudoinverse))
#endif
        call add_test(tests,new_unittest("c_square_pseudoinverse",test_c_square_pseudoinverse))
        call add_test(tests,new_unittest("c_tall_pseudoinverse",test_c_tall_pseudoinverse))
        call add_test(tests,new_unittest("c_wide_pseudoinverse",test_c_wide_pseudoinverse))
        call add_test(tests,new_unittest("c_singular_pseudoinverse",test_c_singular_pseudoinverse))
        call add_test(tests,new_unittest("z_square_pseudoinverse",test_z_square_pseudoinverse))
        call add_test(tests,new_unittest("z_tall_pseudoinverse",test_z_tall_pseudoinverse))
        call add_test(tests,new_unittest("z_wide_pseudoinverse",test_z_wide_pseudoinverse))
        call add_test(tests,new_unittest("z_singular_pseudoinverse",test_z_singular_pseudoinverse))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("y_square_pseudoinverse",test_y_square_pseudoinverse))
        call add_test(tests,new_unittest("y_tall_pseudoinverse",test_y_tall_pseudoinverse))
        call add_test(tests,new_unittest("y_wide_pseudoinverse",test_y_wide_pseudoinverse))
        call add_test(tests,new_unittest("y_singular_pseudoinverse",test_y_singular_pseudoinverse))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("w_square_pseudoinverse",test_w_square_pseudoinverse))
        call add_test(tests,new_unittest("w_tall_pseudoinverse",test_w_tall_pseudoinverse))
        call add_test(tests,new_unittest("w_wide_pseudoinverse",test_w_wide_pseudoinverse))
        call add_test(tests,new_unittest("w_singular_pseudoinverse",test_w_singular_pseudoinverse))
#endif

    end subroutine test_pseudoinverse_matrix

    !> Invert identity matrix
    subroutine test_s_eye_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))

        real(sp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_sp,0.0_sp,i == j)
        end do

        !> Invert funrtion
        inva = pinv(a,err=state)
        
        call check(error,state%ok(),'s pseudoinverse (eye, function): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'s pseudoinverse (eye, function): data convergence')
        if (allocated(error)) return
        
        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        
        call check(error,state%ok(),'s pseudoinverse (eye, subroutine): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'s pseudoinverse (eye, subroutine): data convergence')
        if (allocated(error)) return
        
        !> Operator
        inva = .pinv.a
        
        call check(error,all(abs(a - inva) < tol),'s pseudoinverse (eye, operator): data convergence')
        if (allocated(error)) return

    end subroutine test_s_eye_pseudoinverse

    subroutine test_d_eye_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))

        real(dp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_dp,0.0_dp,i == j)
        end do

        !> Invert funrtion
        inva = pinv(a,err=state)
        
        call check(error,state%ok(),'d pseudoinverse (eye, function): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'d pseudoinverse (eye, function): data convergence')
        if (allocated(error)) return
        
        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        
        call check(error,state%ok(),'d pseudoinverse (eye, subroutine): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'d pseudoinverse (eye, subroutine): data convergence')
        if (allocated(error)) return
        
        !> Operator
        inva = .pinv.a
        
        call check(error,all(abs(a - inva) < tol),'d pseudoinverse (eye, operator): data convergence')
        if (allocated(error)) return

    end subroutine test_d_eye_pseudoinverse

#ifdef LA_WITH_XDP
    subroutine test_x_eye_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))

        real(xdp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_xdp,0.0_xdp,i == j)
        end do

        !> Invert funrtion
        inva = pinv(a,err=state)
        
        call check(error,state%ok(),'x pseudoinverse (eye, function): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'x pseudoinverse (eye, function): data convergence')
        if (allocated(error)) return
        
        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        
        call check(error,state%ok(),'x pseudoinverse (eye, subroutine): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'x pseudoinverse (eye, subroutine): data convergence')
        if (allocated(error)) return
        
        !> Operator
        inva = .pinv.a
        
        call check(error,all(abs(a - inva) < tol),'x pseudoinverse (eye, operator): data convergence')
        if (allocated(error)) return

    end subroutine test_x_eye_pseudoinverse
#endif

#ifdef LA_WITH_QP
    subroutine test_q_eye_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j
        integer(ilp),parameter :: n = 15_ilp
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))

        real(qp) :: a(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge(1.0_qp,0.0_qp,i == j)
        end do

        !> Invert funrtion
        inva = pinv(a,err=state)
        
        call check(error,state%ok(),'q pseudoinverse (eye, function): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'q pseudoinverse (eye, function): data convergence')
        if (allocated(error)) return
        
        !> Inverse subroutine
        call pseudoinvert(a,inva,err=state)
        
        call check(error,state%ok(),'q pseudoinverse (eye, subroutine): '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(a - inva) < tol),'q pseudoinverse (eye, subroutine): data convergence')
        if (allocated(error)) return
        
        !> Operator
        inva = .pinv.a
        
        call check(error,all(abs(a - inva) < tol),'q pseudoinverse (eye, operator): data convergence')
        if (allocated(error)) return

    end subroutine test_q_eye_pseudoinverse
#endif

    !> Test edge case: square matrix
    subroutine test_s_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'s pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'s pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'s pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_s_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_s_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'s pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'s pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'s pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_s_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_s_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'s pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'s pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'s pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_s_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_s_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'s pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'s pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'s pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_s_singular_pseudoinverse

    !> Test edge case: square matrix
    subroutine test_d_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'d pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'d pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'d pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_d_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_d_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        real(dp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'d pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'d pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'d pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_d_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_d_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        real(dp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'d pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'d pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'d pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_d_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_d_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'d pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'d pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'d pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_d_singular_pseudoinverse

#ifdef LA_WITH_XDP
    !> Test edge case: square matrix
    subroutine test_x_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'x pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'x pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'x pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_x_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_x_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'x pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'x pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'x pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_x_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_x_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'x pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'x pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'x pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_x_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_x_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'x pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'x pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'x pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_x_singular_pseudoinverse
#endif

#ifdef LA_WITH_QP
    !> Test edge case: square matrix
    subroutine test_q_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'q pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'q pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'q pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_q_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_q_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'q pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'q pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'q pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_q_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_q_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),inva(n,m)
        
        call random_number(a)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'q pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'q pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'q pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_q_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_q_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),inva(n,n)
        
        call random_number(a)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'q pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'q pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'q pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_q_singular_pseudoinverse
#endif

    !> Test edge case: square matrix
    subroutine test_c_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),inva(n,n)
        real(sp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=sp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'c pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'c pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'c pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_c_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_c_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),inva(n,m)
        real(sp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=sp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'c pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'c pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'c pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_c_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_c_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),inva(n,m)
        real(sp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=sp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'c pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'c pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'c pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_c_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_c_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(sp),parameter :: tol = 1000*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),inva(n,n)
        real(sp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=sp)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'c pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'c pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'c pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_c_singular_pseudoinverse

    !> Test edge case: square matrix
    subroutine test_z_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),inva(n,n)
        real(dp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=dp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'z pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'z pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'z pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_z_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_z_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(m,n),inva(n,m)
        real(dp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=dp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'z pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'z pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'z pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_z_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_z_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(m,n),inva(n,m)
        real(dp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=dp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'z pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'z pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'z pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_z_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_z_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(dp),parameter :: tol = 1000*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),inva(n,n)
        real(dp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=dp)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'z pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'z pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'z pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_z_singular_pseudoinverse

#ifdef LA_WITH_XDP
    !> Test edge case: square matrix
    subroutine test_y_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),inva(n,n)
        real(xdp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=xdp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'y pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'y pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'y pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_y_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_y_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(m,n),inva(n,m)
        real(xdp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=xdp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'y pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'y pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'y pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_y_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_y_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(m,n),inva(n,m)
        real(xdp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=xdp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'y pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'y pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'y pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_y_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_y_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(xdp),parameter :: tol = 1000*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),inva(n,n)
        real(xdp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=xdp)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'y pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'y pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'y pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_y_singular_pseudoinverse
#endif

#ifdef LA_WITH_QP
    !> Test edge case: square matrix
    subroutine test_w_square_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),inva(n,n)
        real(qp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=qp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'w pseudoinverse (square): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'w pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'w pseudoinverse (square, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_w_square_pseudoinverse

    !> Test edge case: tall matrix
    subroutine test_w_tall_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 20,n = 10
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),inva(n,m)
        real(qp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=qp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'w pseudoinverse (tall): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'w pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'w pseudoinverse (tall, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_w_tall_pseudoinverse

    !> Test edge case: wide matrix
    subroutine test_w_wide_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: m = 10,n = 20
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),inva(n,m)
        real(qp) :: rea(m,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=qp)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'w pseudoinverse (wide): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'w pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'w pseudoinverse (wide, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_w_wide_pseudoinverse

    !> Test edge case: singular matrix
    subroutine test_w_singular_pseudoinverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: failed
        integer(ilp),parameter :: n = 10
        real(qp),parameter :: tol = 1000*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),inva(n,n)
        real(qp) :: rea(n,n,2)
        
        call random_number(rea)
        a = cmplx(rea(:,:,1),rea(:,:,2),kind=qp)
        
        ! Make the matrix singular
        a(:,1) = a(:,2)
        
        inva = pinv(a,err=state)
        call check(error,state%ok(),'w pseudoinverse (singular): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(a - matmul(a,matmul(inva,a))) > tol)
        call check(error,failed == 0,'w pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return
        
        failed = count(abs(inva - matmul(inva,matmul(a,inva))) > tol)
        call check(error,failed == 0,'w pseudoinverse (singular, convergence): '//state%print())
        if (allocated(error)) return

    end subroutine test_w_singular_pseudoinverse
#endif

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type),allocatable,intent(inout) :: tests(:)
        type(unittest_type),intent(in) :: new_test
        
        integer :: n
        type(unittest_type),allocatable :: new_tests(:)
        
        if (allocated(tests)) then
            n = size(tests)
        else
            n = 0
        end if
        
        allocate (new_tests(n + 1))
        if (n > 0) new_tests(1:n) = tests(1:n)
                 new_tests(1 + n) = new_test
        call move_alloc(from=new_tests,to=tests)
        
    end subroutine add_test

end module test_linalg_pseudoinverse

program test_pinv
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_pseudoinverse,only:test_pseudoinverse_matrix
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_pseudoinverse",test_pseudoinverse_matrix) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_pinv

