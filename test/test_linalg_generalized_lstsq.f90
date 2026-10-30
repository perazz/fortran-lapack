! Test the generalized least squares solver
module test_linalg_generalized_least_squares
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:generalized_lstsq,solve_generalized_lstsq,generalized_lstsq_space, &
                              lstsq,operator(.inv.)
    use la_state_type,only:la_state,LINALG_ERROR,LINALG_VALUE_ERROR

    implicit none(type,external)
    private

    public :: test_generalized_least_squares

    integer(ilp),parameter :: m = 4,n = 2

    contains

    !> Solve sample generalized least squares problems
    subroutine test_generalized_least_squares(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("generalized_lstsq_identity_s",test_generalized_lstsq_identity_s))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_s",test_generalized_lstsq_correlated_s))
        call add_test(tests,new_unittest("solve_generalized_lstsq_s",test_solve_generalized_lstsq_s))
        call add_test(tests,new_unittest("generalized_lstsq_errors_s",test_generalized_lstsq_errors_s))
        call add_test(tests,new_unittest("generalized_lstsq_identity_d",test_generalized_lstsq_identity_d))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_d",test_generalized_lstsq_correlated_d))
        call add_test(tests,new_unittest("solve_generalized_lstsq_d",test_solve_generalized_lstsq_d))
        call add_test(tests,new_unittest("generalized_lstsq_errors_d",test_generalized_lstsq_errors_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("generalized_lstsq_identity_x",test_generalized_lstsq_identity_x))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_x",test_generalized_lstsq_correlated_x))
        call add_test(tests,new_unittest("solve_generalized_lstsq_x",test_solve_generalized_lstsq_x))
        call add_test(tests,new_unittest("generalized_lstsq_errors_x",test_generalized_lstsq_errors_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("generalized_lstsq_identity_q",test_generalized_lstsq_identity_q))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_q",test_generalized_lstsq_correlated_q))
        call add_test(tests,new_unittest("solve_generalized_lstsq_q",test_solve_generalized_lstsq_q))
        call add_test(tests,new_unittest("generalized_lstsq_errors_q",test_generalized_lstsq_errors_q))
#endif
        call add_test(tests,new_unittest("generalized_lstsq_identity_c",test_generalized_lstsq_identity_c))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_c",test_generalized_lstsq_correlated_c))
        call add_test(tests,new_unittest("solve_generalized_lstsq_c",test_solve_generalized_lstsq_c))
        call add_test(tests,new_unittest("generalized_lstsq_errors_c",test_generalized_lstsq_errors_c))
        call add_test(tests,new_unittest("generalized_lstsq_identity_z",test_generalized_lstsq_identity_z))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_z",test_generalized_lstsq_correlated_z))
        call add_test(tests,new_unittest("solve_generalized_lstsq_z",test_solve_generalized_lstsq_z))
        call add_test(tests,new_unittest("generalized_lstsq_errors_z",test_generalized_lstsq_errors_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("generalized_lstsq_identity_y",test_generalized_lstsq_identity_y))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_y",test_generalized_lstsq_correlated_y))
        call add_test(tests,new_unittest("solve_generalized_lstsq_y",test_solve_generalized_lstsq_y))
        call add_test(tests,new_unittest("generalized_lstsq_errors_y",test_generalized_lstsq_errors_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("generalized_lstsq_identity_w",test_generalized_lstsq_identity_w))
        call add_test(tests,new_unittest("generalized_lstsq_correlated_w",test_generalized_lstsq_correlated_w))
        call add_test(tests,new_unittest("solve_generalized_lstsq_w",test_solve_generalized_lstsq_w))
        call add_test(tests,new_unittest("generalized_lstsq_errors_w",test_generalized_lstsq_errors_w))
#endif

    end subroutine test_generalized_least_squares

    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_s(A,b,L,W,x_ref)
        real(sp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        real(sp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [1.0_sp,2.0_sp,3.0_sp,4.0_sp]
        b = [2.0_sp,4.0_sp,5.0_sp,4.0_sp]
        L(1:4,1) = [2.0_sp,1.0_sp,0.5_sp,0.25_sp]
        L(2:4,2) = [1.5_sp,0.25_sp,-0.5_sp]
        L(3:4,3) = [1.0_sp,0.5_sp]
        L(4,4) = 1.25_sp
        W = matmul(L,transpose(L))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_s

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(sp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_s(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_s

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        type(la_state) :: state
        real(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(sp),allocatable :: x(:)

        call problem_s(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_s

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        type(la_state) :: state
        real(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        real(sp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_s(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_s

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(sp),allocatable :: x(:)

        call problem_s(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_s

    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_d(A,b,L,W,x_ref)
        real(dp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        real(dp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        b = [2.0_dp,4.0_dp,5.0_dp,4.0_dp]
        L(1:4,1) = [2.0_dp,1.0_dp,0.5_dp,0.25_dp]
        L(2:4,2) = [1.5_dp,0.25_dp,-0.5_dp]
        L(3:4,3) = [1.0_dp,0.5_dp]
        L(4,4) = 1.25_dp
        W = matmul(L,transpose(L))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_d

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(dp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_d(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_d

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        type(la_state) :: state
        real(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(dp),allocatable :: x(:)

        call problem_d(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_d

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        type(la_state) :: state
        real(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        real(dp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_d(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_d

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(dp),allocatable :: x(:)

        call problem_d(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_d

#ifdef LA_WITH_XDP
    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_x(A,b,L,W,x_ref)
        real(xdp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        real(xdp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [1.0_xdp,2.0_xdp,3.0_xdp,4.0_xdp]
        b = [2.0_xdp,4.0_xdp,5.0_xdp,4.0_xdp]
        L(1:4,1) = [2.0_xdp,1.0_xdp,0.5_xdp,0.25_xdp]
        L(2:4,2) = [1.5_xdp,0.25_xdp,-0.5_xdp]
        L(3:4,3) = [1.0_xdp,0.5_xdp]
        L(4,4) = 1.25_xdp
        W = matmul(L,transpose(L))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_x

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        real(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(xdp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_x(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_x

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        type(la_state) :: state
        real(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(xdp),allocatable :: x(:)

        call problem_x(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_x

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        type(la_state) :: state
        real(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        real(xdp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_x(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_x

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_x(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(xdp),allocatable :: x(:)

        call problem_x(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_x
#endif

#ifdef LA_WITH_QP
    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_q(A,b,L,W,x_ref)
        real(qp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        real(qp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [1.0_qp,2.0_qp,3.0_qp,4.0_qp]
        b = [2.0_qp,4.0_qp,5.0_qp,4.0_qp]
        L(1:4,1) = [2.0_qp,1.0_qp,0.5_qp,0.25_qp]
        L(2:4,2) = [1.5_qp,0.25_qp,-0.5_qp]
        L(3:4,3) = [1.0_qp,0.5_qp]
        L(4,4) = 1.25_qp
        W = matmul(L,transpose(L))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_q

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(qp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_q(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_q

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        type(la_state) :: state
        real(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(qp),allocatable :: x(:)

        call problem_q(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_q

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        type(la_state) :: state
        real(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        real(qp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_q(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_q

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        real(qp),allocatable :: x(:)

        call problem_q(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_q
#endif

    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_c(A,b,L,W,x_ref)
        complex(sp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        complex(sp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [(1.0_sp,1.0_sp), (2.0_sp,-1.0_sp), (3.0_sp,0.5_sp), (4.0_sp,2.0_sp)]
        b = [(2.0_sp,1.0_sp), (4.0_sp,0.0_sp), (5.0_sp,-1.0_sp), (4.0_sp,3.0_sp)]
        L(1:4,1) = [(2.0_sp,0.0_sp), (1.0_sp,0.5_sp), (0.5_sp,0.0_sp), (0.25_sp,0.5_sp)]
        L(2:4,2) = [(1.5_sp,0.0_sp), (0.25_sp,-0.75_sp), (-0.5_sp,0.0_sp)]
        L(3:4,3) = [(1.0_sp,0.0_sp), (0.5_sp,0.25_sp)]
        L(4,4) = 1.25_sp
        W = matmul(L,conjg(transpose(L)))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_c

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        complex(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(sp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_c(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_c

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        type(la_state) :: state
        complex(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(sp),allocatable :: x(:)

        call problem_c(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_c

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        type(la_state) :: state
        complex(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        complex(sp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_c(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_c

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_c(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        complex(sp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(sp),allocatable :: x(:)

        call problem_c(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_c

    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_z(A,b,L,W,x_ref)
        complex(dp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        complex(dp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [(1.0_dp,1.0_dp), (2.0_dp,-1.0_dp), (3.0_dp,0.5_dp), (4.0_dp,2.0_dp)]
        b = [(2.0_dp,1.0_dp), (4.0_dp,0.0_dp), (5.0_dp,-1.0_dp), (4.0_dp,3.0_dp)]
        L(1:4,1) = [(2.0_dp,0.0_dp), (1.0_dp,0.5_dp), (0.5_dp,0.0_dp), (0.25_dp,0.5_dp)]
        L(2:4,2) = [(1.5_dp,0.0_dp), (0.25_dp,-0.75_dp), (-0.5_dp,0.0_dp)]
        L(3:4,3) = [(1.0_dp,0.0_dp), (0.5_dp,0.25_dp)]
        L(4,4) = 1.25_dp
        W = matmul(L,conjg(transpose(L)))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_z

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        complex(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(dp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_z(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_z

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        type(la_state) :: state
        complex(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(dp),allocatable :: x(:)

        call problem_z(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_z

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        type(la_state) :: state
        complex(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        complex(dp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_z(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_z

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_z(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        complex(dp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(dp),allocatable :: x(:)

        call problem_z(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_z

#ifdef LA_WITH_XDP
    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_y(A,b,L,W,x_ref)
        complex(xdp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        complex(xdp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [(1.0_xdp,1.0_xdp), (2.0_xdp,-1.0_xdp), (3.0_xdp,0.5_xdp), (4.0_xdp,2.0_xdp)]
        b = [(2.0_xdp,1.0_xdp), (4.0_xdp,0.0_xdp), (5.0_xdp,-1.0_xdp), (4.0_xdp,3.0_xdp)]
        L(1:4,1) = [(2.0_xdp,0.0_xdp), (1.0_xdp,0.5_xdp), (0.5_xdp,0.0_xdp), (0.25_xdp,0.5_xdp)]
        L(2:4,2) = [(1.5_xdp,0.0_xdp), (0.25_xdp,-0.75_xdp), (-0.5_xdp,0.0_xdp)]
        L(3:4,3) = [(1.0_xdp,0.0_xdp), (0.5_xdp,0.25_xdp)]
        L(4,4) = 1.25_xdp
        W = matmul(L,conjg(transpose(L)))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_y

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        complex(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(xdp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_y(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_y

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        type(la_state) :: state
        complex(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(xdp),allocatable :: x(:)

        call problem_y(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_y

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        type(la_state) :: state
        complex(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        complex(xdp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_y(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_y

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_y(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        complex(xdp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(xdp),allocatable :: x(:)

        call problem_y(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_y
#endif

#ifdef LA_WITH_QP
    !> A straight-line fit, the lower triangular square root L of its covariance W = L*L^H,
    !> and the solution of the whitened problem min ||L^{-1}(b - A*x)||
    subroutine problem_w(A,b,L,W,x_ref)
        complex(qp),intent(out) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)

        complex(qp) :: Linv(m,m),LA(m,n),Lb(m)

        A(:,1) = 1
        L = 0
        A(:,2) = [(1.0_qp,1.0_qp), (2.0_qp,-1.0_qp), (3.0_qp,0.5_qp), (4.0_qp,2.0_qp)]
        b = [(2.0_qp,1.0_qp), (4.0_qp,0.0_qp), (5.0_qp,-1.0_qp), (4.0_qp,3.0_qp)]
        L(1:4,1) = [(2.0_qp,0.0_qp), (1.0_qp,0.5_qp), (0.5_qp,0.0_qp), (0.25_qp,0.5_qp)]
        L(2:4,2) = [(1.5_qp,0.0_qp), (0.25_qp,-0.75_qp), (-0.5_qp,0.0_qp)]
        L(3:4,3) = [(1.0_qp,0.0_qp), (0.5_qp,0.25_qp)]
        L(4,4) = 1.25_qp
        W = matmul(L,conjg(transpose(L)))

        Linv = .inv.L
        LA = matmul(Linv,A)
        Lb = matmul(Linv,b)
        x_ref = lstsq(LA,Lb)

    end subroutine problem_w

    !> With an identity covariance, the generalized solution is the ordinary one
    subroutine test_generalized_lstsq_identity_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        complex(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(qp),allocatable :: x(:)
        integer(ilp) :: i

        call problem_w(A,b,L,W,x_ref)

        W = 0
        do i = 1,m
           W(i,i) = 1
        end do

        x = generalized_lstsq(W,A,b)

        call check(error,all(abs(x - lstsq(A,b)) < tol),'generalized_lstsq with W = I is not lstsq')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_identity_w

    !> Correlated errors, from the covariance or from any of its square roots
    subroutine test_generalized_lstsq_correlated_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        type(la_state) :: state
        complex(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(qp),allocatable :: x(:)

        call problem_w(A,b,L,W,x_ref)

        x = generalized_lstsq(W,A,b)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the covariance')
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq with err')
        if (allocated(error)) return

        x = generalized_lstsq(L,A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from the Cholesky factor')
        if (allocated(error)) return

        ! L*P, with P a permutation, is a square root of W that is not triangular
        x = generalized_lstsq(L(:, [3,1,4,2]),A,b,prefactored_w=.true.)
        call check(error,all(abs(x - x_ref) < tol),'generalized_lstsq from a non-triangular square root')
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_correlated_w

    !> Subroutine form: inputs kept by default, caller workspace, in-place solve
    subroutine test_solve_generalized_lstsq_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        type(la_state) :: state
        complex(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n),x(n),A0(m,n),W0(m,m)
        complex(qp),allocatable :: work(:)
        integer(ilp) :: lwork

        call problem_w(A,b,L,W,x_ref)
        A0 = A
        W0 = W

        call solve_generalized_lstsq(W,A,b,x,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return
        call check(error,all(A == A0) .and. all(W == W0),'solve_generalized_lstsq changed a or w')
        if (allocated(error)) return

        call generalized_lstsq_space(A,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        allocate (work(lwork))
        call solve_generalized_lstsq(W,A,b,x,storage=work,overwrite_a=.true.,overwrite_w=.true.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        call check(error,all(abs(x - x_ref) < tol),'solve_generalized_lstsq in place')
        if (allocated(error)) return

        A = A0
        W = W0
        call solve_generalized_lstsq(W,A,b,x,storage=work(:lwork - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'short workspace returned '//state%print())
        if (allocated(error)) return

    end subroutine test_solve_generalized_lstsq_w

    !> Invalid and singular problems
    subroutine test_generalized_lstsq_errors_w(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        complex(qp) :: A(m,n),b(m),L(m,m),W(m,m),x_ref(n)
        complex(qp),allocatable :: x(:)

        call problem_w(A,b,L,W,x_ref)

        x = generalized_lstsq(-W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'negative definite w returned '//state%print())
        if (allocated(error)) return
        call check(error,size(x) == 0,'a failed solve returns an empty solution')
        if (allocated(error)) return

        x = generalized_lstsq(W(1:n,1:n),transpose(A),b(1:n),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'a with m < n returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W(1:m - 1,1:m - 1),A,b,err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'w of the wrong size returned '//state%print())
        if (allocated(error)) return

        x = generalized_lstsq(W,A,b(1:m - 1),err=state)
        call check(error,state%state == LINALG_VALUE_ERROR,'b of the wrong size returned '//state%print())
        if (allocated(error)) return

        A(:,2) = 0
        x = generalized_lstsq(W,A,b,err=state)
        call check(error,state%state == LINALG_ERROR,'rank-deficient a returned '//state%print())
        if (allocated(error)) return

    end subroutine test_generalized_lstsq_errors_w
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

end module test_linalg_generalized_least_squares

program test_generalized_lstsq
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_generalized_least_squares,only:test_generalized_least_squares
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_generalized_least_squares",test_generalized_least_squares) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_generalized_lstsq
