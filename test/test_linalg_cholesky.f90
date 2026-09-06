! Test Cholesky factorization
module test_linalg_cholesky
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:cholesky,chol
    use la_state_type,only:la_state

    implicit none(type,external)
    private
    
    public :: test_cholesky_factorization

    contains

    !> Cholesky factorization tests
    subroutine test_cholesky_factorization(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))
        
        call add_test(tests,new_unittest("least_cholesky_s",test_cholesky_s))
        call add_test(tests,new_unittest("least_cholesky_d",test_cholesky_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("least_cholesky_x",test_cholesky_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("least_cholesky_q",test_cholesky_q))
#endif
        call add_test(tests,new_unittest("least_cholesky_c",test_cholesky_c))
        call add_test(tests,new_unittest("least_cholesky_z",test_cholesky_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("least_cholesky_y",test_cholesky_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("least_cholesky_w",test_cholesky_w))
#endif

    end subroutine test_cholesky_factorization

    !> Cholesky factorization of a random matrix
    subroutine test_cholesky_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_sp,0.0000_sp,0.0000_sp]
        l(2,:) = [6.1237_sp,4.1833_sp,0.0000_sp]
        l(3,:) = [22.4537_sp,20.9165_sp,6.1101_sp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_s

    subroutine test_cholesky_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_dp,0.0000_dp,0.0000_dp]
        l(2,:) = [6.1237_dp,4.1833_dp,0.0000_dp]
        l(3,:) = [22.4537_dp,20.9165_dp,6.1101_dp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_d

#ifdef LA_WITH_XDP
    subroutine test_cholesky_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_xdp,0.0000_xdp,0.0000_xdp]
        l(2,:) = [6.1237_xdp,4.1833_xdp,0.0000_xdp]
        l(3,:) = [22.4537_xdp,20.9165_xdp,6.1101_xdp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_x
#endif

#ifdef LA_WITH_QP
    subroutine test_cholesky_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_qp,0.0000_qp,0.0000_qp]
        l(2,:) = [6.1237_qp,4.1833_qp,0.0000_qp]
        l(3,:) = [22.4537_qp,20.9165_qp,6.1101_qp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_q
#endif

    subroutine test_cholesky_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_sp,0.0000_sp,0.0000_sp]
        l(2,:) = [6.1237_sp,4.1833_sp,0.0000_sp]
        l(3,:) = [22.4537_sp,20.9165_sp,6.1101_sp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_c

    subroutine test_cholesky_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_dp,0.0000_dp,0.0000_dp]
        l(2,:) = [6.1237_dp,4.1833_dp,0.0000_dp]
        l(3,:) = [22.4537_dp,20.9165_dp,6.1101_dp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_z

#ifdef LA_WITH_XDP
    subroutine test_cholesky_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_xdp,0.0000_xdp,0.0000_xdp]
        l(2,:) = [6.1237_xdp,4.1833_xdp,0.0000_xdp]
        l(3,:) = [22.4537_xdp,20.9165_xdp,6.1101_xdp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_y
#endif

#ifdef LA_WITH_QP
    subroutine test_cholesky_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),l(n,n)
        type(la_state) :: state
        
        ! Set real matrix
        a(1,:) = [6,15,55]
        a(2,:) = [15,55,225]
        a(3,:) = [55,225,979]
        
        ! Set result (lower factor)
        l(1,:) = [2.4495_qp,0.0000_qp,0.0000_qp]
        l(2,:) = [6.1237_qp,4.1833_qp,0.0000_qp]
        l(3,:) = [22.4537_qp,20.9165_qp,6.1101_qp]
        
        ! 1) Cholesky factorization with full matrices
        call cholesky(a,l,other_zeroed=.true.,err=state)
        
        call check(error,state%ok(),'cholesky (subr) :: '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (subr) :: data converged')
        if (allocated(error)) return
        
        ! 2) Function interface
        l = chol(a,other_zeroed=.true.)
        
        call check(error,all(abs(a - matmul(l,transpose(l))) < tol),'cholesky (function) :: data converged')
        if (allocated(error)) return
        
    end subroutine test_cholesky_w
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

end module test_linalg_cholesky

program test_cholesky
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_cholesky,only:test_cholesky_factorization
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_cholesky",test_cholesky_factorization) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_cholesky
