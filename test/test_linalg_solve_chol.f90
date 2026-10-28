! Test solve_chol, solve_lower_chol, solve_upper_chol
module test_linalg_solve_chol
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:cholesky,solve_chol,solve_lower_chol,solve_upper_chol
    use la_state_type,only:la_state

    implicit none(type,external)
    private
    
    public :: test_solve_chol_factorization

    contains

    !> Cholesky solve tests
    subroutine test_solve_chol_factorization(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))
        
        call add_test(tests,new_unittest("solve_lower_chol_s",test_solve_lower_chol_s))
        call add_test(tests,new_unittest("solve_upper_chol_s",test_solve_upper_chol_s))
        call add_test(tests,new_unittest("solve_chol_s",test_solve_chol_s))
        call add_test(tests,new_unittest("solve_chol_upper_s",test_solve_chol_upper_s))
        call add_test(tests,new_unittest("solve_chol_overwrite_s",test_solve_chol_overwrite_s))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_s",test_solve_lower_chol_multi_rhs_s))
        call add_test(tests,new_unittest("solve_chol_indefinite_s",test_solve_chol_indefinite_s))
        call add_test(tests,new_unittest("solve_chol_semidefinite_s",test_solve_chol_semidefinite_s))
        call add_test(tests,new_unittest("solve_lower_chol_d",test_solve_lower_chol_d))
        call add_test(tests,new_unittest("solve_upper_chol_d",test_solve_upper_chol_d))
        call add_test(tests,new_unittest("solve_chol_d",test_solve_chol_d))
        call add_test(tests,new_unittest("solve_chol_upper_d",test_solve_chol_upper_d))
        call add_test(tests,new_unittest("solve_chol_overwrite_d",test_solve_chol_overwrite_d))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_d",test_solve_lower_chol_multi_rhs_d))
        call add_test(tests,new_unittest("solve_chol_indefinite_d",test_solve_chol_indefinite_d))
        call add_test(tests,new_unittest("solve_chol_semidefinite_d",test_solve_chol_semidefinite_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("solve_lower_chol_x",test_solve_lower_chol_x))
        call add_test(tests,new_unittest("solve_upper_chol_x",test_solve_upper_chol_x))
        call add_test(tests,new_unittest("solve_chol_x",test_solve_chol_x))
        call add_test(tests,new_unittest("solve_chol_upper_x",test_solve_chol_upper_x))
        call add_test(tests,new_unittest("solve_chol_overwrite_x",test_solve_chol_overwrite_x))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_x",test_solve_lower_chol_multi_rhs_x))
        call add_test(tests,new_unittest("solve_chol_indefinite_x",test_solve_chol_indefinite_x))
        call add_test(tests,new_unittest("solve_chol_semidefinite_x",test_solve_chol_semidefinite_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("solve_lower_chol_q",test_solve_lower_chol_q))
        call add_test(tests,new_unittest("solve_upper_chol_q",test_solve_upper_chol_q))
        call add_test(tests,new_unittest("solve_chol_q",test_solve_chol_q))
        call add_test(tests,new_unittest("solve_chol_upper_q",test_solve_chol_upper_q))
        call add_test(tests,new_unittest("solve_chol_overwrite_q",test_solve_chol_overwrite_q))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_q",test_solve_lower_chol_multi_rhs_q))
        call add_test(tests,new_unittest("solve_chol_indefinite_q",test_solve_chol_indefinite_q))
        call add_test(tests,new_unittest("solve_chol_semidefinite_q",test_solve_chol_semidefinite_q))
#endif
        call add_test(tests,new_unittest("solve_lower_chol_c",test_solve_lower_chol_c))
        call add_test(tests,new_unittest("solve_upper_chol_c",test_solve_upper_chol_c))
        call add_test(tests,new_unittest("solve_chol_c",test_solve_chol_c))
        call add_test(tests,new_unittest("solve_chol_upper_c",test_solve_chol_upper_c))
        call add_test(tests,new_unittest("solve_chol_overwrite_c",test_solve_chol_overwrite_c))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_c",test_solve_lower_chol_multi_rhs_c))
        call add_test(tests,new_unittest("solve_chol_indefinite_c",test_solve_chol_indefinite_c))
        call add_test(tests,new_unittest("solve_chol_semidefinite_c",test_solve_chol_semidefinite_c))
        call add_test(tests,new_unittest("solve_lower_chol_z",test_solve_lower_chol_z))
        call add_test(tests,new_unittest("solve_upper_chol_z",test_solve_upper_chol_z))
        call add_test(tests,new_unittest("solve_chol_z",test_solve_chol_z))
        call add_test(tests,new_unittest("solve_chol_upper_z",test_solve_chol_upper_z))
        call add_test(tests,new_unittest("solve_chol_overwrite_z",test_solve_chol_overwrite_z))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_z",test_solve_lower_chol_multi_rhs_z))
        call add_test(tests,new_unittest("solve_chol_indefinite_z",test_solve_chol_indefinite_z))
        call add_test(tests,new_unittest("solve_chol_semidefinite_z",test_solve_chol_semidefinite_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("solve_lower_chol_y",test_solve_lower_chol_y))
        call add_test(tests,new_unittest("solve_upper_chol_y",test_solve_upper_chol_y))
        call add_test(tests,new_unittest("solve_chol_y",test_solve_chol_y))
        call add_test(tests,new_unittest("solve_chol_upper_y",test_solve_chol_upper_y))
        call add_test(tests,new_unittest("solve_chol_overwrite_y",test_solve_chol_overwrite_y))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_y",test_solve_lower_chol_multi_rhs_y))
        call add_test(tests,new_unittest("solve_chol_indefinite_y",test_solve_chol_indefinite_y))
        call add_test(tests,new_unittest("solve_chol_semidefinite_y",test_solve_chol_semidefinite_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("solve_lower_chol_w",test_solve_lower_chol_w))
        call add_test(tests,new_unittest("solve_upper_chol_w",test_solve_upper_chol_w))
        call add_test(tests,new_unittest("solve_chol_w",test_solve_chol_w))
        call add_test(tests,new_unittest("solve_chol_upper_w",test_solve_chol_upper_w))
        call add_test(tests,new_unittest("solve_chol_overwrite_w",test_solve_chol_overwrite_w))
        call add_test(tests,new_unittest("solve_lower_chol_multi_rhs_w",test_solve_lower_chol_multi_rhs_w))
        call add_test(tests,new_unittest("solve_chol_indefinite_w",test_solve_chol_indefinite_w))
        call add_test(tests,new_unittest("solve_chol_semidefinite_w",test_solve_chol_semidefinite_w))
#endif

    end subroutine test_solve_chol_factorization

    !> Test solve_lower_chol with pre-computed lower Cholesky factors
    subroutine test_solve_lower_chol_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_s(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_s
    
    subroutine test_solve_lower_chol_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_d(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_lower_chol_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_x(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_lower_chol_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_q(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_q
#endif
    
    subroutine test_solve_lower_chol_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_c(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_c
    
    subroutine test_solve_lower_chol_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_z(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_lower_chol_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_y(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_lower_chol_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),l(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol: solution mismatch: '// &
                   mismatch_w(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_w
#endif
    
    !> Test solve_upper_chol with pre-computed upper Cholesky factors
    subroutine test_solve_upper_chol_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_s(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_s
    
    subroutine test_solve_upper_chol_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_d(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_upper_chol_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_x(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_upper_chol_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_q(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_q
#endif
    
    subroutine test_solve_upper_chol_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_c(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_c
    
    subroutine test_solve_upper_chol_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_z(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_upper_chol_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_y(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_upper_chol_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),u(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! Compute Cholesky factorization (upper triangular: A = U^T * U)
        call cholesky(a,u,lower=.false.,err=state)
        
        call check(error,state%ok(),'cholesky factorization (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using upper Cholesky factors
        call solve_upper_chol(u,b,x,err=state)
        
        call check(error,state%ok(),'solve_upper_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_upper_chol: solution mismatch: '// &
                   mismatch_w(a,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]), &
                                   tol,state,factor=u,lower=.false.))
        if (allocated(error)) return
        
    end subroutine test_solve_upper_chol_w
#endif
    
    !> Test solve_chol (one-shot) - default preserves A
    subroutine test_solve_chol_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_s(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_s
    
    subroutine test_solve_chol_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_d(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_x(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_q(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_q
#endif
    
    subroutine test_solve_chol_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_c(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_c
    
    subroutine test_solve_chol_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_z(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_y(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with default overwrite_a=.false. (A should be preserved)
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol: solution mismatch: '// &
                   mismatch_w(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol: A should be preserved by default')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_w
#endif
    
    !> Test solve_chol with upper triangular (lower=.false.)
    subroutine test_solve_chol_upper_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_s(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_s
    
    subroutine test_solve_chol_upper_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_d(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_upper_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_x(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_upper_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_q(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_q
#endif
    
    subroutine test_solve_chol_upper_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_c(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_c
    
    subroutine test_solve_chol_upper_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_z(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_upper_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_y(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_upper_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with upper triangular (A should be preserved by default)
        call solve_chol(a,b,x,lower=.false.,err=state)
        
        call check(error,state%ok(),'solve_chol (upper) failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol (upper): solution mismatch: '// &
                   mismatch_w(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A is preserved (default behavior)
        call check(error,all(abs(a - a_copy) < tol),'solve_chol (upper): A should be preserved')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_upper_w
#endif
    
    !> Test solve_chol with overwrite_a=.true.
    subroutine test_solve_chol_overwrite_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_s(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_s
    
    subroutine test_solve_chol_overwrite_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_d(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_overwrite_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_x(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_overwrite_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_q(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_q
#endif
    
    subroutine test_solve_chol_overwrite_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_c(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_c
    
    subroutine test_solve_chol_overwrite_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_z(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_overwrite_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_y(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_overwrite_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),a_copy(n,n),b(n),x(n),x_expected(n)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        a_copy = a
        
        ! Known solution
        x_expected = [1,2,3]
        
        ! Compute RHS: b = A * x_expected
        b = matmul(a,x_expected)
        
        ! One-shot solve with overwrite_a=.true. (A will be destroyed)
        call solve_chol(a,b,x,lower=.true.,overwrite_a=.true.,err=state)
        
        call check(error,state%ok(),'solve_chol overwrite failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_chol overwrite: solution mismatch: '// &
                   mismatch_w(a_copy,reshape(b, [n,1_ilp]),reshape(x, [n,1_ilp]),reshape(x_expected, [n,1_ilp]),tol,state))
        if (allocated(error)) return
        
        ! Check that A was overwritten (not equal to original)
        call check(error,any(abs(a - a_copy) > tol),'solve_chol: A should be overwritten with overwrite_a=.true.')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_overwrite_w
#endif
    
    !> Test solve_lower_chol with multiple RHS
    subroutine test_solve_lower_chol_multi_rhs_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_s(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_s
    
    subroutine test_solve_lower_chol_multi_rhs_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_d(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_lower_chol_multi_rhs_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_x(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_lower_chol_multi_rhs_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_q(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_q
#endif
    
    subroutine test_solve_lower_chol_multi_rhs_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_c(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_c
    
    subroutine test_solve_lower_chol_multi_rhs_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_z(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_lower_chol_multi_rhs_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_y(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_lower_chol_multi_rhs_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        integer(ilp),parameter :: nrhs = 2_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),l(n,n),b(n,nrhs),x(n,nrhs),x_expected(n,nrhs)
        type(la_state) :: state
        
        ! Set symmetric positive definite matrix
        a(1,:) = [4,2,2]
        a(2,:) = [2,5,1]
        a(3,:) = [2,1,6]
        
        ! Known solutions (two RHS)
        x_expected(:,1) = [1,2,3]
        x_expected(:,2) = [4,5,6]
        
        ! Compute RHS: B = A * X_expected. gfortran 15.1 and 15.2, and 16 snapshots before 16.1,
        ! miscompile this complex matmul at -O3 with FMA (GCC PR 124802), so they get it written out
#if defined(__GFORTRAN__) && ((__GNUC__ == 15 && __GNUC_MINOR__ < 3) || (__GNUC__ == 16 && __GNUC_MINOR__ == 0))
        b(:,1) = [14,15,22]
        b(:,2) = [38,39,49]
#else
        b = matmul(a,x_expected)
#endif
        
        ! Compute Cholesky factorization
        call cholesky(a,l,lower=.true.,err=state)
        
        call check(error,state%ok(),'cholesky factorization failed: '//state%print())
        if (allocated(error)) return
        
        ! Solve using lower Cholesky factors
        call solve_lower_chol(l,b,x,err=state)
        
        call check(error,state%ok(),'solve_lower_chol multi-rhs failed: '//state%print())
        if (allocated(error)) return
        
        ! Check solution
        call check(error,all(abs(x - x_expected) < tol),'solve_lower_chol multi-rhs: solution mismatch: '// &
                   mismatch_w(a,b,x,x_expected,tol,state,factor=l,lower=.true.))
        if (allocated(error)) return
        
    end subroutine test_solve_lower_chol_multi_rhs_w
#endif
    
    !> Test solve_chol with symmetric indefinite matrix (should fail)
    subroutine test_solve_chol_indefinite_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(sp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_s
    
    subroutine test_solve_chol_indefinite_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(dp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_indefinite_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(xdp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_indefinite_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(qp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_q
#endif
    
    subroutine test_solve_chol_indefinite_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(sp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_c
    
    subroutine test_solve_chol_indefinite_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(dp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_indefinite_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(xdp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_indefinite_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(qp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric INDEFINITE matrix (eigenvalues: 3, -1)
        a(1,:) = [1,2]
        a(2,:) = [2,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for indefinite matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not positive definite)
        call check(error,state%error(),'solve_chol should fail for indefinite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_indefinite_w
#endif
    
    !> Test solve_chol with symmetric positive semi-definite matrix (should fail)
    subroutine test_solve_chol_semidefinite_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(sp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_s
    
    subroutine test_solve_chol_semidefinite_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(dp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_d
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_semidefinite_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(xdp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_x
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_semidefinite_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        real(qp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_q
#endif
    
    subroutine test_solve_chol_semidefinite_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(sp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_c
    
    subroutine test_solve_chol_semidefinite_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(dp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_z
    
#ifdef LA_WITH_XDP
    subroutine test_solve_chol_semidefinite_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(xdp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_y
#endif
    
#ifdef LA_WITH_QP
    subroutine test_solve_chol_semidefinite_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 2_ilp
        complex(qp) :: a(n,n),b(n),x(n)
        type(la_state) :: state
        
        ! Set symmetric positive SEMI-DEFINITE matrix (eigenvalues: 2, 0 - singular)
        a(1,:) = [1,1]
        a(2,:) = [1,1]
        
        ! Arbitrary RHS
        b = [1,1]
        
        ! solve_chol should fail for semi-definite (singular) matrix
        call solve_chol(a,b,x,lower=.true.,err=state)
        
        ! Check that it failed (not strictly positive definite)
        call check(error,state%error(),'solve_chol should fail for positive semi-definite matrix')
        if (allocated(error)) return
        
    end subroutine test_solve_chol_semidefinite_w
#endif
    
    !> What a failed solution check saw: the largest error, where it sits and the tolerance, the
    !> residual of the system, the state the solver returned, the computed and expected
    !> solutions and, when the test factorized first, how far the factor is from the matrix
    function mismatch_s(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        real(sp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        real(sp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(sp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        real(sp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*(g0,:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        real(sp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,transpose(factor))
            else
                rebuilt = matmul(transpose(factor),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_s

    function mismatch_d(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        real(dp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        real(dp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(dp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        real(dp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*(g0,:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        real(dp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,transpose(factor))
            else
                rebuilt = matmul(transpose(factor),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_d

#ifdef LA_WITH_XDP
    function mismatch_x(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        real(xdp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        real(xdp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(xdp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        real(xdp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*(g0,:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        real(xdp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,transpose(factor))
            else
                rebuilt = matmul(transpose(factor),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_x
#endif

#ifdef LA_WITH_QP
    function mismatch_q(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        real(qp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        real(qp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(qp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        real(qp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*(g0,:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        real(qp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,transpose(factor))
            else
                rebuilt = matmul(transpose(factor),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_q
#endif

    function mismatch_c(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        complex(sp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        complex(sp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(sp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        complex(sp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*("(",g0,",",g0,")",:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        complex(sp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,conjg(transpose(factor)))
            else
                rebuilt = matmul(conjg(transpose(factor)),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_c

    function mismatch_z(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        complex(dp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        complex(dp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(dp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        complex(dp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*("(",g0,",",g0,")",:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        complex(dp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,conjg(transpose(factor)))
            else
                rebuilt = matmul(conjg(transpose(factor)),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_z

#ifdef LA_WITH_XDP
    function mismatch_y(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        complex(xdp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        complex(xdp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(xdp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        complex(xdp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*("(",g0,",",g0,")",:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        complex(xdp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,conjg(transpose(factor)))
            else
                rebuilt = matmul(conjg(transpose(factor)),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_y
#endif

#ifdef LA_WITH_QP
    function mismatch_w(a,b,x,x_expected,tol,state,factor,lower) result(msg)
        !> Matrix of the system, as it was before the solve
        complex(qp),intent(in) :: a(:,:)
        !> Right-hand sides, computed and expected solutions, one column each
        complex(qp),intent(in) :: b(:,:),x(:,:),x_expected(:,:)
        !> Tolerance the check held each entry of the solution to
        real(qp),intent(in) :: tol
        !> State the solver returned
        type(la_state),intent(in) :: state
        !> [optional] Cholesky factor the solve used, lower or upper as `lower` says (default lower)
        complex(qp),optional,intent(in) :: factor(:,:)
        logical,optional,intent(in) :: lower
        character(:),allocatable :: msg

        character(*),parameter :: values = '(*("(",g0,",",g0,")",:,1x))'
        character(len=2048) :: buffer
        integer :: loc(2)
        logical :: lower_
        complex(qp),allocatable :: rebuilt(:,:)

        loc = maxloc(abs(x - x_expected))
        write (buffer,'(a,es10.3,a,i0,a,i0,a,es10.3,a,es10.3)') 'max|x-x_expected|=', &
            maxval(abs(x - x_expected)),' at (',loc(1),',',loc(2),'), tol=',tol, &
            ', max|A*x-b|=',maxval(abs(matmul(a,x) - b))
        msg = trim(buffer)

        if (present(factor)) then
            lower_ = .true.
            if (present(lower)) lower_ = lower
            if (lower_) then
                rebuilt = matmul(factor,conjg(transpose(factor)))
            else
                rebuilt = matmul(conjg(transpose(factor)),factor)
            end if
            write (buffer,'(a,es10.3)') ', max|factor product - A|=',maxval(abs(rebuilt - a))
            msg = msg//trim(buffer)
        end if

        msg = msg//', solver state: '//state%print()
        write (buffer,values) x
        msg = msg//', x=['//trim(buffer)//']'
        write (buffer,values) x_expected
        msg = msg//', x_expected=['//trim(buffer)//']'

    end function mismatch_w
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

end module test_linalg_solve_chol

program test_solve_chol
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_solve_chol,only:test_solve_chol_factorization
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_solve_chol",test_solve_chol_factorization) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_solve_chol
