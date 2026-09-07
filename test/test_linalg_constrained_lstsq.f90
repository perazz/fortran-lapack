! Test least squares solver
module test_linalg_constrained_least_squares
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:constrained_lstsq,solve_constrained_lstsq,constrained_lstsq_space
    use la_state_type,only:la_state

    implicit none(type,external)
    private
    
    public :: test_constrained_least_squares

    contains

    !> Solve sample least squares problems
    subroutine test_constrained_least_squares(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)
        
        allocate (tests(0))
        
        call add_test(tests,new_unittest("constrained_least_squares_randm_s",test_constrained_lstsq_random_s))
        call add_test(tests,new_unittest("constrained_least_squares_randm_d",test_constrained_lstsq_random_d))
        call add_test(tests,new_unittest("constrained_least_squares_randm_q",test_constrained_lstsq_random_q))
        call add_test(tests,new_unittest("constrained_least_squares_randm_c",test_constrained_lstsq_random_c))
        call add_test(tests,new_unittest("constrained_least_squares_randm_z",test_constrained_lstsq_random_z))
        call add_test(tests,new_unittest("constrained_least_squares_randm_w",test_constrained_lstsq_random_w))

    end subroutine test_constrained_least_squares
    
    !> Fit from random array
    subroutine test_constrained_lstsq_random_s(error)
        type(error_type),allocatable,intent(out) :: error
        type(la_state) :: state
        integer(ilp),parameter :: m = 5,n = 4,p = 3
        !> Least-squares cost.
        real(sp) :: A(m,n),b(m)
        !> Equality constraints.
        real(sp) :: C(p,n),d(p)
        !> Solution.
        real(sp) :: x(n),x_true(n)
        !> Workspace.
        integer(ilp) :: lwork
        real(sp),allocatable :: work(:)

        !> Least-squares cost.
        A(1,:) = [1.0_sp,1.0_sp,1.0_sp,1.0_sp]
        A(2,:) = [1.0_sp,3.0_sp,1.0_sp,1.0_sp]
        A(3,:) = [1.0_sp,-1.0_sp,3.0_sp,1.0_sp]
        A(4,:) = [1.0_sp,1.0_sp,1.0_sp,3.0_sp]
        A(5,:) = [1.0_sp,1.0_sp,1.0_sp,-1.0_sp]

        b = [2.0_sp,1.0_sp,6.0_sp,3.0_sp,1.0_sp]

        !> Equality constraints.
        C(1,:) = [1.0_sp,1.0_sp,1.0_sp,-1.0_sp]
        C(2,:) = [1.0_sp,-1.0_sp,1.0_sp,1.0_sp]
        C(3,:) = [1.0_sp,1.0_sp,-1.0_sp,1.0_sp]

        d = [1.0_sp,3.0_sp,-1.0_sp]

        !----- Function interface -----
        x = constrained_lstsq(A,b,C,d,err=state)
        x_true = [0.5_sp,-0.5_sp,1.5_sp,0.5_sp]

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

        !----- Subroutine interface -----
        call solve_constrained_lstsq(A,b,C,d,x,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

        !----- Pre-allocated storage -----
        call constrained_lstsq_space(A,C,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        allocate (work(lwork))
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

        !----- Overwrite matrices (performance) -----
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,overwrite_matrices=.true.,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

    end subroutine test_constrained_lstsq_random_s
    
    !> Fit from random array
    subroutine test_constrained_lstsq_random_d(error)
        type(error_type),allocatable,intent(out) :: error
        type(la_state) :: state
        integer(ilp),parameter :: m = 5,n = 4,p = 3
        !> Least-squares cost.
        real(dp) :: A(m,n),b(m)
        !> Equality constraints.
        real(dp) :: C(p,n),d(p)
        !> Solution.
        real(dp) :: x(n),x_true(n)
        !> Workspace.
        integer(ilp) :: lwork
        real(dp),allocatable :: work(:)

        !> Least-squares cost.
        A(1,:) = [1.0_dp,1.0_dp,1.0_dp,1.0_dp]
        A(2,:) = [1.0_dp,3.0_dp,1.0_dp,1.0_dp]
        A(3,:) = [1.0_dp,-1.0_dp,3.0_dp,1.0_dp]
        A(4,:) = [1.0_dp,1.0_dp,1.0_dp,3.0_dp]
        A(5,:) = [1.0_dp,1.0_dp,1.0_dp,-1.0_dp]

        b = [2.0_dp,1.0_dp,6.0_dp,3.0_dp,1.0_dp]

        !> Equality constraints.
        C(1,:) = [1.0_dp,1.0_dp,1.0_dp,-1.0_dp]
        C(2,:) = [1.0_dp,-1.0_dp,1.0_dp,1.0_dp]
        C(3,:) = [1.0_dp,1.0_dp,-1.0_dp,1.0_dp]

        d = [1.0_dp,3.0_dp,-1.0_dp]

        !----- Function interface -----
        x = constrained_lstsq(A,b,C,d,err=state)
        x_true = [0.5_dp,-0.5_dp,1.5_dp,0.5_dp]

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

        !----- Subroutine interface -----
        call solve_constrained_lstsq(A,b,C,d,x,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

        !----- Pre-allocated storage -----
        call constrained_lstsq_space(A,C,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        allocate (work(lwork))
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

        !----- Overwrite matrices (performance) -----
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,overwrite_matrices=.true.,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

    end subroutine test_constrained_lstsq_random_d
    
    !> Fit from random array
    subroutine test_constrained_lstsq_random_q(error)
        type(error_type),allocatable,intent(out) :: error
        type(la_state) :: state
        integer(ilp),parameter :: m = 5,n = 4,p = 3
        !> Least-squares cost.
        real(qp) :: A(m,n),b(m)
        !> Equality constraints.
        real(qp) :: C(p,n),d(p)
        !> Solution.
        real(qp) :: x(n),x_true(n)
        !> Workspace.
        integer(ilp) :: lwork
        real(qp),allocatable :: work(:)

        !> Least-squares cost.
        A(1,:) = [1.0_qp,1.0_qp,1.0_qp,1.0_qp]
        A(2,:) = [1.0_qp,3.0_qp,1.0_qp,1.0_qp]
        A(3,:) = [1.0_qp,-1.0_qp,3.0_qp,1.0_qp]
        A(4,:) = [1.0_qp,1.0_qp,1.0_qp,3.0_qp]
        A(5,:) = [1.0_qp,1.0_qp,1.0_qp,-1.0_qp]

        b = [2.0_qp,1.0_qp,6.0_qp,3.0_qp,1.0_qp]

        !> Equality constraints.
        C(1,:) = [1.0_qp,1.0_qp,1.0_qp,-1.0_qp]
        C(2,:) = [1.0_qp,-1.0_qp,1.0_qp,1.0_qp]
        C(3,:) = [1.0_qp,1.0_qp,-1.0_qp,1.0_qp]

        d = [1.0_qp,3.0_qp,-1.0_qp]

        !----- Function interface -----
        x = constrained_lstsq(A,b,C,d,err=state)
        x_true = [0.5_qp,-0.5_qp,1.5_qp,0.5_qp]

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

        !----- Subroutine interface -----
        call solve_constrained_lstsq(A,b,C,d,x,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

        !----- Pre-allocated storage -----
        call constrained_lstsq_space(A,C,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        allocate (work(lwork))
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

        !----- Overwrite matrices (performance) -----
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,overwrite_matrices=.true.,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

    end subroutine test_constrained_lstsq_random_q
    
    !> Fit from random array
    subroutine test_constrained_lstsq_random_c(error)
        type(error_type),allocatable,intent(out) :: error
        type(la_state) :: state
        integer(ilp),parameter :: m = 5,n = 4,p = 3
        !> Least-squares cost.
        complex(sp) :: A(m,n),b(m)
        !> Equality constraints.
        complex(sp) :: C(p,n),d(p)
        !> Solution.
        complex(sp) :: x(n),x_true(n)
        !> Workspace.
        integer(ilp) :: lwork
        complex(sp),allocatable :: work(:)

        !> Least-squares cost.
        A(1,:) = [1.0_sp,1.0_sp,1.0_sp,1.0_sp]
        A(2,:) = [1.0_sp,3.0_sp,1.0_sp,1.0_sp]
        A(3,:) = [1.0_sp,-1.0_sp,3.0_sp,1.0_sp]
        A(4,:) = [1.0_sp,1.0_sp,1.0_sp,3.0_sp]
        A(5,:) = [1.0_sp,1.0_sp,1.0_sp,-1.0_sp]

        b = [2.0_sp,1.0_sp,6.0_sp,3.0_sp,1.0_sp]

        !> Equality constraints.
        C(1,:) = [1.0_sp,1.0_sp,1.0_sp,-1.0_sp]
        C(2,:) = [1.0_sp,-1.0_sp,1.0_sp,1.0_sp]
        C(3,:) = [1.0_sp,1.0_sp,-1.0_sp,1.0_sp]

        d = [1.0_sp,3.0_sp,-1.0_sp]

        !----- Function interface -----
        x = constrained_lstsq(A,b,C,d,err=state)
        x_true = [0.5_sp,-0.5_sp,1.5_sp,0.5_sp]

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

        !----- Subroutine interface -----
        call solve_constrained_lstsq(A,b,C,d,x,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

        !----- Pre-allocated storage -----
        call constrained_lstsq_space(A,C,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        allocate (work(lwork))
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

        !----- Overwrite matrices (performance) -----
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,overwrite_matrices=.true.,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_sp),'Solver converged')
        if (allocated(error)) return

    end subroutine test_constrained_lstsq_random_c
    
    !> Fit from random array
    subroutine test_constrained_lstsq_random_z(error)
        type(error_type),allocatable,intent(out) :: error
        type(la_state) :: state
        integer(ilp),parameter :: m = 5,n = 4,p = 3
        !> Least-squares cost.
        complex(dp) :: A(m,n),b(m)
        !> Equality constraints.
        complex(dp) :: C(p,n),d(p)
        !> Solution.
        complex(dp) :: x(n),x_true(n)
        !> Workspace.
        integer(ilp) :: lwork
        complex(dp),allocatable :: work(:)

        !> Least-squares cost.
        A(1,:) = [1.0_dp,1.0_dp,1.0_dp,1.0_dp]
        A(2,:) = [1.0_dp,3.0_dp,1.0_dp,1.0_dp]
        A(3,:) = [1.0_dp,-1.0_dp,3.0_dp,1.0_dp]
        A(4,:) = [1.0_dp,1.0_dp,1.0_dp,3.0_dp]
        A(5,:) = [1.0_dp,1.0_dp,1.0_dp,-1.0_dp]

        b = [2.0_dp,1.0_dp,6.0_dp,3.0_dp,1.0_dp]

        !> Equality constraints.
        C(1,:) = [1.0_dp,1.0_dp,1.0_dp,-1.0_dp]
        C(2,:) = [1.0_dp,-1.0_dp,1.0_dp,1.0_dp]
        C(3,:) = [1.0_dp,1.0_dp,-1.0_dp,1.0_dp]

        d = [1.0_dp,3.0_dp,-1.0_dp]

        !----- Function interface -----
        x = constrained_lstsq(A,b,C,d,err=state)
        x_true = [0.5_dp,-0.5_dp,1.5_dp,0.5_dp]

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

        !----- Subroutine interface -----
        call solve_constrained_lstsq(A,b,C,d,x,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

        !----- Pre-allocated storage -----
        call constrained_lstsq_space(A,C,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        allocate (work(lwork))
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

        !----- Overwrite matrices (performance) -----
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,overwrite_matrices=.true.,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_dp),'Solver converged')
        if (allocated(error)) return

    end subroutine test_constrained_lstsq_random_z
    
    !> Fit from random array
    subroutine test_constrained_lstsq_random_w(error)
        type(error_type),allocatable,intent(out) :: error
        type(la_state) :: state
        integer(ilp),parameter :: m = 5,n = 4,p = 3
        !> Least-squares cost.
        complex(qp) :: A(m,n),b(m)
        !> Equality constraints.
        complex(qp) :: C(p,n),d(p)
        !> Solution.
        complex(qp) :: x(n),x_true(n)
        !> Workspace.
        integer(ilp) :: lwork
        complex(qp),allocatable :: work(:)

        !> Least-squares cost.
        A(1,:) = [1.0_qp,1.0_qp,1.0_qp,1.0_qp]
        A(2,:) = [1.0_qp,3.0_qp,1.0_qp,1.0_qp]
        A(3,:) = [1.0_qp,-1.0_qp,3.0_qp,1.0_qp]
        A(4,:) = [1.0_qp,1.0_qp,1.0_qp,3.0_qp]
        A(5,:) = [1.0_qp,1.0_qp,1.0_qp,-1.0_qp]

        b = [2.0_qp,1.0_qp,6.0_qp,3.0_qp,1.0_qp]

        !> Equality constraints.
        C(1,:) = [1.0_qp,1.0_qp,1.0_qp,-1.0_qp]
        C(2,:) = [1.0_qp,-1.0_qp,1.0_qp,1.0_qp]
        C(3,:) = [1.0_qp,1.0_qp,-1.0_qp,1.0_qp]

        d = [1.0_qp,3.0_qp,-1.0_qp]

        !----- Function interface -----
        x = constrained_lstsq(A,b,C,d,err=state)
        x_true = [0.5_qp,-0.5_qp,1.5_qp,0.5_qp]

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

        !----- Subroutine interface -----
        call solve_constrained_lstsq(A,b,C,d,x,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

        !----- Pre-allocated storage -----
        call constrained_lstsq_space(A,C,lwork,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        allocate (work(lwork))
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

        !----- Overwrite matrices (performance) -----
        call solve_constrained_lstsq(A,b,C,d,x,storage=work,overwrite_matrices=.true.,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < 1.0e-4_qp),'Solver converged')
        if (allocated(error)) return

    end subroutine test_constrained_lstsq_random_w
    
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

end module test_linalg_constrained_least_squares

program test_constrained_lstsq
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_constrained_least_squares,only:test_constrained_least_squares
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_constrained_least_squares",test_constrained_least_squares) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_constrained_lstsq
