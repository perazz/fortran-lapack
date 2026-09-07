! Test least squares solver
module test_linalg_least_squares
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:lstsq,lstsq_space,solve_lstsq,weighted_lstsq,solve_weighted_lstsq
    use la_state_type,only:la_state

    implicit none(type,external)
    private
    
    public :: test_least_squares

    contains

    !> Solve sample least squares problems
    subroutine test_least_squares(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)
        
        allocate (tests(0))
        
        call add_test(tests,new_unittest("issue_823",test_issue_823))

        call add_test(tests,new_unittest("least_squares_s",test_lstsq_one_s))
        call add_test(tests,new_unittest("least_squares_randm_s",test_lstsq_random_s))
        call add_test(tests,new_unittest("weighted_lstsq_s",test_weighted_lstsq_s))
        call add_test(tests,new_unittest("weighted_lstsq_effect_s",test_weighted_lstsq_effect_s))
        call add_test(tests,new_unittest("weighted_lstsq_negative_s",test_weighted_lstsq_negative_s))
        call add_test(tests,new_unittest("solve_weighted_lstsq_s",test_solve_weighted_lstsq_s))
        call add_test(tests,new_unittest("solve_weighted_lstsq_effect_s",test_solve_weighted_lstsq_effect_s))
        call add_test(tests,new_unittest("solve_weighted_lstsq_negative_s",test_solve_weighted_lstsq_negative_s))
        call add_test(tests,new_unittest("least_squares_d",test_lstsq_one_d))
        call add_test(tests,new_unittest("least_squares_randm_d",test_lstsq_random_d))
        call add_test(tests,new_unittest("weighted_lstsq_d",test_weighted_lstsq_d))
        call add_test(tests,new_unittest("weighted_lstsq_effect_d",test_weighted_lstsq_effect_d))
        call add_test(tests,new_unittest("weighted_lstsq_negative_d",test_weighted_lstsq_negative_d))
        call add_test(tests,new_unittest("solve_weighted_lstsq_d",test_solve_weighted_lstsq_d))
        call add_test(tests,new_unittest("solve_weighted_lstsq_effect_d",test_solve_weighted_lstsq_effect_d))
        call add_test(tests,new_unittest("solve_weighted_lstsq_negative_d",test_solve_weighted_lstsq_negative_d))
        call add_test(tests,new_unittest("least_squares_q",test_lstsq_one_q))
        call add_test(tests,new_unittest("least_squares_randm_q",test_lstsq_random_q))
        call add_test(tests,new_unittest("weighted_lstsq_q",test_weighted_lstsq_q))
        call add_test(tests,new_unittest("weighted_lstsq_effect_q",test_weighted_lstsq_effect_q))
        call add_test(tests,new_unittest("weighted_lstsq_negative_q",test_weighted_lstsq_negative_q))
        call add_test(tests,new_unittest("solve_weighted_lstsq_q",test_solve_weighted_lstsq_q))
        call add_test(tests,new_unittest("solve_weighted_lstsq_effect_q",test_solve_weighted_lstsq_effect_q))
        call add_test(tests,new_unittest("solve_weighted_lstsq_negative_q",test_solve_weighted_lstsq_negative_q))

        call add_test(tests,new_unittest("weighted_lstsq_c",test_weighted_lstsq_c))
        call add_test(tests,new_unittest("solve_weighted_lstsq_c",test_solve_weighted_lstsq_c))
        call add_test(tests,new_unittest("weighted_lstsq_z",test_weighted_lstsq_z))
        call add_test(tests,new_unittest("solve_weighted_lstsq_z",test_solve_weighted_lstsq_z))
        call add_test(tests,new_unittest("weighted_lstsq_w",test_weighted_lstsq_w))
        call add_test(tests,new_unittest("solve_weighted_lstsq_w",test_solve_weighted_lstsq_w))

        call add_test(tests,new_unittest("lstsq_space_s",test_lstsq_space_s))
        call add_test(tests,new_unittest("lstsq_space_d",test_lstsq_space_d))
        call add_test(tests,new_unittest("lstsq_space_q",test_lstsq_space_q))
        call add_test(tests,new_unittest("lstsq_space_c",test_lstsq_space_c))
        call add_test(tests,new_unittest("lstsq_space_z",test_lstsq_space_z))
        call add_test(tests,new_unittest("lstsq_space_w",test_lstsq_space_w))

    end subroutine test_least_squares

    !> Solve into a caller-provided solution array and caller-provided working arrays
    subroutine test_lstsq_space_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 6,n = 3
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))

        integer(ilp) :: i,j,rank,lrwork,liwork
        integer(ilp),allocatable,target :: iwork(:)
        real(sp),allocatable,target :: rwork(:)
        real(sp) :: singvals(min(m,n))
        real(sp) :: a(m,n),b(m),x(n),x_true(n)
        real(sp),allocatable :: x_alloc(:)

        ! Vandermonde system with an exact solution
        do j = 1,n
            x_true(j) = j
            do i = 1,m
                a(i,j) = i**(j - 1)
            end do
        end do
        b = matmul(a,x_true)

        call lstsq_space(a,b,lrwork,liwork)
        allocate (rwork(lrwork),iwork(liwork))
        call solve_lstsq(a,b,x,real_storage=rwork,int_storage=iwork, &
                         singvals=singvals,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < tol),'solve_lstsq into caller storage')
        if (allocated(error)) return

        call check(error,rank == n,'solve_lstsq returns the full rank')
        if (allocated(error)) return

        call check(error,all(singvals(1:n - 1) >= singvals(2:n)),'singular values in decreasing order')
        if (allocated(error)) return

        ! Working space that is too small is a value error
        call solve_lstsq(a,b,x,real_storage=rwork(1:1),int_storage=iwork,err=state)

        call check(error,state%error(),'insufficient working space is an error')
        if (allocated(error)) return

        ! A solution array that is too small is a value error
        call solve_lstsq(a,b,x(1:n - 1),err=state)

        call check(error,state%error(),'a solution array that is too short is an error')
        if (allocated(error)) return

        ! The function form returns an empty solution for an invalid problem
        x_alloc = lstsq(a,b(1:m - 1),err=state)

        call check(error,state%error() .and. size(x_alloc) == 0,'an invalid problem returns an empty solution')
        if (allocated(error)) return

    end subroutine test_lstsq_space_s
    !> Solve into a caller-provided solution array and caller-provided working arrays
    subroutine test_lstsq_space_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 6,n = 3
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))

        integer(ilp) :: i,j,rank,lrwork,liwork
        integer(ilp),allocatable,target :: iwork(:)
        real(dp),allocatable,target :: rwork(:)
        real(dp) :: singvals(min(m,n))
        real(dp) :: a(m,n),b(m),x(n),x_true(n)
        real(dp),allocatable :: x_alloc(:)

        ! Vandermonde system with an exact solution
        do j = 1,n
            x_true(j) = j
            do i = 1,m
                a(i,j) = i**(j - 1)
            end do
        end do
        b = matmul(a,x_true)

        call lstsq_space(a,b,lrwork,liwork)
        allocate (rwork(lrwork),iwork(liwork))
        call solve_lstsq(a,b,x,real_storage=rwork,int_storage=iwork, &
                         singvals=singvals,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < tol),'solve_lstsq into caller storage')
        if (allocated(error)) return

        call check(error,rank == n,'solve_lstsq returns the full rank')
        if (allocated(error)) return

        call check(error,all(singvals(1:n - 1) >= singvals(2:n)),'singular values in decreasing order')
        if (allocated(error)) return

        ! Working space that is too small is a value error
        call solve_lstsq(a,b,x,real_storage=rwork(1:1),int_storage=iwork,err=state)

        call check(error,state%error(),'insufficient working space is an error')
        if (allocated(error)) return

        ! A solution array that is too small is a value error
        call solve_lstsq(a,b,x(1:n - 1),err=state)

        call check(error,state%error(),'a solution array that is too short is an error')
        if (allocated(error)) return

        ! The function form returns an empty solution for an invalid problem
        x_alloc = lstsq(a,b(1:m - 1),err=state)

        call check(error,state%error() .and. size(x_alloc) == 0,'an invalid problem returns an empty solution')
        if (allocated(error)) return

    end subroutine test_lstsq_space_d
    !> Solve into a caller-provided solution array and caller-provided working arrays
    subroutine test_lstsq_space_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 6,n = 3
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))

        integer(ilp) :: i,j,rank,lrwork,liwork
        integer(ilp),allocatable,target :: iwork(:)
        real(qp),allocatable,target :: rwork(:)
        real(qp) :: singvals(min(m,n))
        real(qp) :: a(m,n),b(m),x(n),x_true(n)
        real(qp),allocatable :: x_alloc(:)

        ! Vandermonde system with an exact solution
        do j = 1,n
            x_true(j) = j
            do i = 1,m
                a(i,j) = i**(j - 1)
            end do
        end do
        b = matmul(a,x_true)

        call lstsq_space(a,b,lrwork,liwork)
        allocate (rwork(lrwork),iwork(liwork))
        call solve_lstsq(a,b,x,real_storage=rwork,int_storage=iwork, &
                         singvals=singvals,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < tol),'solve_lstsq into caller storage')
        if (allocated(error)) return

        call check(error,rank == n,'solve_lstsq returns the full rank')
        if (allocated(error)) return

        call check(error,all(singvals(1:n - 1) >= singvals(2:n)),'singular values in decreasing order')
        if (allocated(error)) return

        ! Working space that is too small is a value error
        call solve_lstsq(a,b,x,real_storage=rwork(1:1),int_storage=iwork,err=state)

        call check(error,state%error(),'insufficient working space is an error')
        if (allocated(error)) return

        ! A solution array that is too small is a value error
        call solve_lstsq(a,b,x(1:n - 1),err=state)

        call check(error,state%error(),'a solution array that is too short is an error')
        if (allocated(error)) return

        ! The function form returns an empty solution for an invalid problem
        x_alloc = lstsq(a,b(1:m - 1),err=state)

        call check(error,state%error() .and. size(x_alloc) == 0,'an invalid problem returns an empty solution')
        if (allocated(error)) return

    end subroutine test_lstsq_space_q
    !> Solve into a caller-provided solution array and caller-provided working arrays
    subroutine test_lstsq_space_c(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 6,n = 3
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))

        integer(ilp) :: i,j,rank,lrwork,liwork
        integer(ilp),allocatable,target :: iwork(:)
        real(sp),allocatable,target :: rwork(:)
        real(sp) :: singvals(min(m,n))
        complex(sp) :: a(m,n),b(m),x(n),x_true(n)
        complex(sp),allocatable :: x_alloc(:)
        integer(ilp) :: lcwork
        complex(sp),allocatable,target :: cwork(:)

        ! Vandermonde system with an exact solution
        do j = 1,n
            x_true(j) = j
            do i = 1,m
                a(i,j) = i**(j - 1)
            end do
        end do
        b = matmul(a,x_true)

        call lstsq_space(a,b,lrwork,liwork,lcwork)
        allocate (rwork(lrwork),iwork(liwork),cwork(lcwork))
        call solve_lstsq(a,b,x,real_storage=rwork,int_storage=iwork,cmpl_storage=cwork, &
                         singvals=singvals,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < tol),'solve_lstsq into caller storage')
        if (allocated(error)) return

        call check(error,rank == n,'solve_lstsq returns the full rank')
        if (allocated(error)) return

        call check(error,all(singvals(1:n - 1) >= singvals(2:n)),'singular values in decreasing order')
        if (allocated(error)) return

        ! Working space that is too small is a value error
        call solve_lstsq(a,b,x,real_storage=rwork(1:1),int_storage=iwork,cmpl_storage=cwork,err=state)

        call check(error,state%error(),'insufficient working space is an error')
        if (allocated(error)) return

        ! A solution array that is too small is a value error
        call solve_lstsq(a,b,x(1:n - 1),err=state)

        call check(error,state%error(),'a solution array that is too short is an error')
        if (allocated(error)) return

        ! The function form returns an empty solution for an invalid problem
        x_alloc = lstsq(a,b(1:m - 1),err=state)

        call check(error,state%error() .and. size(x_alloc) == 0,'an invalid problem returns an empty solution')
        if (allocated(error)) return

    end subroutine test_lstsq_space_c
    !> Solve into a caller-provided solution array and caller-provided working arrays
    subroutine test_lstsq_space_z(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 6,n = 3
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))

        integer(ilp) :: i,j,rank,lrwork,liwork
        integer(ilp),allocatable,target :: iwork(:)
        real(dp),allocatable,target :: rwork(:)
        real(dp) :: singvals(min(m,n))
        complex(dp) :: a(m,n),b(m),x(n),x_true(n)
        complex(dp),allocatable :: x_alloc(:)
        integer(ilp) :: lcwork
        complex(dp),allocatable,target :: cwork(:)

        ! Vandermonde system with an exact solution
        do j = 1,n
            x_true(j) = j
            do i = 1,m
                a(i,j) = i**(j - 1)
            end do
        end do
        b = matmul(a,x_true)

        call lstsq_space(a,b,lrwork,liwork,lcwork)
        allocate (rwork(lrwork),iwork(liwork),cwork(lcwork))
        call solve_lstsq(a,b,x,real_storage=rwork,int_storage=iwork,cmpl_storage=cwork, &
                         singvals=singvals,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < tol),'solve_lstsq into caller storage')
        if (allocated(error)) return

        call check(error,rank == n,'solve_lstsq returns the full rank')
        if (allocated(error)) return

        call check(error,all(singvals(1:n - 1) >= singvals(2:n)),'singular values in decreasing order')
        if (allocated(error)) return

        ! Working space that is too small is a value error
        call solve_lstsq(a,b,x,real_storage=rwork(1:1),int_storage=iwork,cmpl_storage=cwork,err=state)

        call check(error,state%error(),'insufficient working space is an error')
        if (allocated(error)) return

        ! A solution array that is too small is a value error
        call solve_lstsq(a,b,x(1:n - 1),err=state)

        call check(error,state%error(),'a solution array that is too short is an error')
        if (allocated(error)) return

        ! The function form returns an empty solution for an invalid problem
        x_alloc = lstsq(a,b(1:m - 1),err=state)

        call check(error,state%error() .and. size(x_alloc) == 0,'an invalid problem returns an empty solution')
        if (allocated(error)) return

    end subroutine test_lstsq_space_z
    !> Solve into a caller-provided solution array and caller-provided working arrays
    subroutine test_lstsq_space_w(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 6,n = 3
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))

        integer(ilp) :: i,j,rank,lrwork,liwork
        integer(ilp),allocatable,target :: iwork(:)
        real(qp),allocatable,target :: rwork(:)
        real(qp) :: singvals(min(m,n))
        complex(qp) :: a(m,n),b(m),x(n),x_true(n)
        complex(qp),allocatable :: x_alloc(:)
        integer(ilp) :: lcwork
        complex(qp),allocatable,target :: cwork(:)

        ! Vandermonde system with an exact solution
        do j = 1,n
            x_true(j) = j
            do i = 1,m
                a(i,j) = i**(j - 1)
            end do
        end do
        b = matmul(a,x_true)

        call lstsq_space(a,b,lrwork,liwork,lcwork)
        allocate (rwork(lrwork),iwork(liwork),cwork(lcwork))
        call solve_lstsq(a,b,x,real_storage=rwork,int_storage=iwork,cmpl_storage=cwork, &
                         singvals=singvals,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,all(abs(x - x_true) < tol),'solve_lstsq into caller storage')
        if (allocated(error)) return

        call check(error,rank == n,'solve_lstsq returns the full rank')
        if (allocated(error)) return

        call check(error,all(singvals(1:n - 1) >= singvals(2:n)),'singular values in decreasing order')
        if (allocated(error)) return

        ! Working space that is too small is a value error
        call solve_lstsq(a,b,x,real_storage=rwork(1:1),int_storage=iwork,cmpl_storage=cwork,err=state)

        call check(error,state%error(),'insufficient working space is an error')
        if (allocated(error)) return

        ! A solution array that is too small is a value error
        call solve_lstsq(a,b,x(1:n - 1),err=state)

        call check(error,state%error(),'a solution array that is too short is an error')
        if (allocated(error)) return

        ! The function form returns an empty solution for an invalid problem
        x_alloc = lstsq(a,b(1:m - 1),err=state)

        call check(error,state%error() .and. size(x_alloc) == 0,'an invalid problem returns an empty solution')
        if (allocated(error)) return

    end subroutine test_lstsq_space_w
    
    !> Simple polynomial fit
    subroutine test_lstsq_one_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp) :: rank

        !> Example scattered data
        real(sp),parameter :: x(*) = real([1.0,2.5,3.5,4.0,5.0,7.0,8.5],sp)
        real(sp),parameter :: y(*) = real([0.3,1.1,1.5,2.0,3.2,6.6,8.6],sp)
        real(sp),parameter :: ab(*) = real([0.20925829,0.12013861],sp)

        real(sp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(p - ab) < 1.0e-4_sp),'data converged')
        if (allocated(error)) return
        
        call check(error,rank == 2,'matrix rank == 2')
        if (allocated(error)) return

    end subroutine test_lstsq_one_s
    
    !> Fit from random array
    subroutine test_lstsq_random_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: n = 12,m = 3
        real :: Arnd(n,m),xrnd(m)
        real(sp),allocatable :: x(:)
        real(sp) :: xsol(m),y(n),A(n,m)

        ! Random coefficient matrix and solution
        call random_number(Arnd)
        call random_number(xrnd)
        
        ! Compute rhs
        A = real(Arnd,sp)
        xsol = real(xrnd,sp)
        y = matmul(A,xsol)

        ! Find polynomial
        x = lstsq(A,y,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        ! Check size
        call check(error,size(x) == m)
        if (allocated(error)) return
        
        call check(error,all(abs(x - xsol) < 1.0e-4_sp),'data converged')
        if (allocated(error)) return
        
    end subroutine test_lstsq_random_s
    
    !> Simple polynomial fit
    subroutine test_lstsq_one_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp) :: rank

        !> Example scattered data
        real(dp),parameter :: x(*) = real([1.0,2.5,3.5,4.0,5.0,7.0,8.5],dp)
        real(dp),parameter :: y(*) = real([0.3,1.1,1.5,2.0,3.2,6.6,8.6],dp)
        real(dp),parameter :: ab(*) = real([0.20925829,0.12013861],dp)

        real(dp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(p - ab) < 1.0e-4_dp),'data converged')
        if (allocated(error)) return
        
        call check(error,rank == 2,'matrix rank == 2')
        if (allocated(error)) return

    end subroutine test_lstsq_one_d
    
    !> Fit from random array
    subroutine test_lstsq_random_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: n = 12,m = 3
        real :: Arnd(n,m),xrnd(m)
        real(dp),allocatable :: x(:)
        real(dp) :: xsol(m),y(n),A(n,m)

        ! Random coefficient matrix and solution
        call random_number(Arnd)
        call random_number(xrnd)
        
        ! Compute rhs
        A = real(Arnd,dp)
        xsol = real(xrnd,dp)
        y = matmul(A,xsol)

        ! Find polynomial
        x = lstsq(A,y,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        ! Check size
        call check(error,size(x) == m)
        if (allocated(error)) return
        
        call check(error,all(abs(x - xsol) < 1.0e-4_dp),'data converged')
        if (allocated(error)) return
        
    end subroutine test_lstsq_random_d
    
    !> Simple polynomial fit
    subroutine test_lstsq_one_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp) :: rank

        !> Example scattered data
        real(qp),parameter :: x(*) = real([1.0,2.5,3.5,4.0,5.0,7.0,8.5],qp)
        real(qp),parameter :: y(*) = real([0.3,1.1,1.5,2.0,3.2,6.6,8.6],qp)
        real(qp),parameter :: ab(*) = real([0.20925829,0.12013861],qp)

        real(qp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,rank=rank,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(p - ab) < 1.0e-4_qp),'data converged')
        if (allocated(error)) return
        
        call check(error,rank == 2,'matrix rank == 2')
        if (allocated(error)) return

    end subroutine test_lstsq_one_q
    
    !> Fit from random array
    subroutine test_lstsq_random_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: n = 12,m = 3
        real :: Arnd(n,m),xrnd(m)
        real(qp),allocatable :: x(:)
        real(qp) :: xsol(m),y(n),A(n,m)

        ! Random coefficient matrix and solution
        call random_number(Arnd)
        call random_number(xrnd)
        
        ! Compute rhs
        A = real(Arnd,qp)
        xsol = real(xrnd,qp)
        y = matmul(A,xsol)

        ! Find polynomial
        x = lstsq(A,y,err=state)

        call check(error,state%ok(),state%print())
        if (allocated(error)) return
        
        ! Check size
        call check(error,size(x) == m)
        if (allocated(error)) return
        
        call check(error,all(abs(x - xsol) < 1.0e-4_qp),'data converged')
        if (allocated(error)) return
        
    end subroutine test_lstsq_random_q
    
    !-------------------------------------------------------------
    !-----     Weighted Least-Squares Tests                  -----
    !-------------------------------------------------------------

    !> Test basic weighted least-squares: uniform weights must match ordinary least-squares
    subroutine test_weighted_lstsq_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n),A_original(m,n),A_copy(m,n),b(m)
        real(sp) :: w(m)
        real(sp),allocatable :: x_weighted(:),x_ols(:)

        ! Simple test case
        A(:,1) = 1.0_sp
        A(:,2) = [1.0_sp,2.0_sp,3.0_sp,4.0_sp]
        b = [2.0_sp,4.0_sp,5.0_sp,4.0_sp]
        w = 1.0_sp ! Uniform weights = OLS
        A_original = A ! Save original to verify A is preserved
        A_copy = A ! Save copy for lstsq comparison

        x_weighted = weighted_lstsq(w,A,b,err=state)

        call check(error,state%ok(),'weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error,size(x_weighted) == n,'weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! Verify A is preserved by default (overwrite_a contract)
        call check(error,all(A == A_original), &
                   'weighted_lstsq must preserve A by default')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares exactly
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_weighted - x_ols) < tol), &
                   'weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_s

    !> Test that non-uniform weights change the solution
    subroutine test_weighted_lstsq_effect_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n),b(m)
        real(sp) :: w_uniform(m),w_nonuniform(m)
        real(sp),allocatable :: x_uniform(:),x_weighted(:)

        ! Setup problem
        A(:,1) = 1.0_sp
        A(:,2) = [1.0_sp,2.0_sp,3.0_sp,4.0_sp]
        b = [1.0_sp,3.0_sp,2.0_sp,5.0_sp]
        
        w_uniform = 1.0_sp
        w_nonuniform = [10.0_sp,1.0_sp,1.0_sp,10.0_sp] ! Weight first and last more

        x_uniform = weighted_lstsq(w_uniform,A,b,err=state)
        call check(error,state%ok(),'uniform weighted_lstsq failed')
        if (allocated(error)) return

        x_weighted = weighted_lstsq(w_nonuniform,A,b,err=state)
        call check(error,state%ok(),'non-uniform weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error,any(abs(x_uniform - x_weighted) > tol), &
                   'weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_effect_s

    !> Test error on negative weights
    subroutine test_weighted_lstsq_negative_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(sp) :: A(3,2),b(3)
        real(sp) :: w(3)
        real(sp),allocatable :: x(:)

        A = 1.0_sp
        b = 1.0_sp
        w = [-1.0_sp,1.0_sp,1.0_sp] ! Invalid: negative weight!

        x = weighted_lstsq(w,A,b,err=state)

        call check(error,state%error(),'weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_negative_s

    !> Test basic weighted least-squares: uniform weights must match ordinary least-squares
    subroutine test_weighted_lstsq_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n),A_original(m,n),A_copy(m,n),b(m)
        real(dp) :: w(m)
        real(dp),allocatable :: x_weighted(:),x_ols(:)

        ! Simple test case
        A(:,1) = 1.0_dp
        A(:,2) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        b = [2.0_dp,4.0_dp,5.0_dp,4.0_dp]
        w = 1.0_dp ! Uniform weights = OLS
        A_original = A ! Save original to verify A is preserved
        A_copy = A ! Save copy for lstsq comparison

        x_weighted = weighted_lstsq(w,A,b,err=state)

        call check(error,state%ok(),'weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error,size(x_weighted) == n,'weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! Verify A is preserved by default (overwrite_a contract)
        call check(error,all(A == A_original), &
                   'weighted_lstsq must preserve A by default')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares exactly
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_weighted - x_ols) < tol), &
                   'weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_d

    !> Test that non-uniform weights change the solution
    subroutine test_weighted_lstsq_effect_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n),b(m)
        real(dp) :: w_uniform(m),w_nonuniform(m)
        real(dp),allocatable :: x_uniform(:),x_weighted(:)

        ! Setup problem
        A(:,1) = 1.0_dp
        A(:,2) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        b = [1.0_dp,3.0_dp,2.0_dp,5.0_dp]
        
        w_uniform = 1.0_dp
        w_nonuniform = [10.0_dp,1.0_dp,1.0_dp,10.0_dp] ! Weight first and last more

        x_uniform = weighted_lstsq(w_uniform,A,b,err=state)
        call check(error,state%ok(),'uniform weighted_lstsq failed')
        if (allocated(error)) return

        x_weighted = weighted_lstsq(w_nonuniform,A,b,err=state)
        call check(error,state%ok(),'non-uniform weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error,any(abs(x_uniform - x_weighted) > tol), &
                   'weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_effect_d

    !> Test error on negative weights
    subroutine test_weighted_lstsq_negative_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(dp) :: A(3,2),b(3)
        real(dp) :: w(3)
        real(dp),allocatable :: x(:)

        A = 1.0_dp
        b = 1.0_dp
        w = [-1.0_dp,1.0_dp,1.0_dp] ! Invalid: negative weight!

        x = weighted_lstsq(w,A,b,err=state)

        call check(error,state%error(),'weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_negative_d

    !> Test basic weighted least-squares: uniform weights must match ordinary least-squares
    subroutine test_weighted_lstsq_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: A(m,n),A_original(m,n),A_copy(m,n),b(m)
        real(qp) :: w(m)
        real(qp),allocatable :: x_weighted(:),x_ols(:)

        ! Simple test case
        A(:,1) = 1.0_qp
        A(:,2) = [1.0_qp,2.0_qp,3.0_qp,4.0_qp]
        b = [2.0_qp,4.0_qp,5.0_qp,4.0_qp]
        w = 1.0_qp ! Uniform weights = OLS
        A_original = A ! Save original to verify A is preserved
        A_copy = A ! Save copy for lstsq comparison

        x_weighted = weighted_lstsq(w,A,b,err=state)

        call check(error,state%ok(),'weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error,size(x_weighted) == n,'weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! Verify A is preserved by default (overwrite_a contract)
        call check(error,all(A == A_original), &
                   'weighted_lstsq must preserve A by default')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares exactly
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_weighted - x_ols) < tol), &
                   'weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_q

    !> Test that non-uniform weights change the solution
    subroutine test_weighted_lstsq_effect_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: A(m,n),b(m)
        real(qp) :: w_uniform(m),w_nonuniform(m)
        real(qp),allocatable :: x_uniform(:),x_weighted(:)

        ! Setup problem
        A(:,1) = 1.0_qp
        A(:,2) = [1.0_qp,2.0_qp,3.0_qp,4.0_qp]
        b = [1.0_qp,3.0_qp,2.0_qp,5.0_qp]
        
        w_uniform = 1.0_qp
        w_nonuniform = [10.0_qp,1.0_qp,1.0_qp,10.0_qp] ! Weight first and last more

        x_uniform = weighted_lstsq(w_uniform,A,b,err=state)
        call check(error,state%ok(),'uniform weighted_lstsq failed')
        if (allocated(error)) return

        x_weighted = weighted_lstsq(w_nonuniform,A,b,err=state)
        call check(error,state%ok(),'non-uniform weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error,any(abs(x_uniform - x_weighted) > tol), &
                   'weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_effect_q

    !> Test error on negative weights
    subroutine test_weighted_lstsq_negative_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(qp) :: A(3,2),b(3)
        real(qp) :: w(3)
        real(qp),allocatable :: x(:)

        A = 1.0_qp
        b = 1.0_qp
        w = [-1.0_qp,1.0_qp,1.0_qp] ! Invalid: negative weight!

        x = weighted_lstsq(w,A,b,err=state)

        call check(error,state%error(),'weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_negative_q

    !> Test weighted least-squares with complex data: uniform weights must match lstsq
    subroutine test_weighted_lstsq_c(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: A(m,n),A_copy(m,n),b(m)
        real(sp) :: w(m)
        complex(sp),allocatable :: x_weighted(:),x_ols(:)

        ! Complex test case
        A(:,1) = cmplx(1.0_sp,0.5_sp,kind=sp)
        A(:,2) = cmplx([1.0_sp,2.0_sp,3.0_sp,4.0_sp], &
                       [0.1_sp,0.2_sp,0.3_sp,0.4_sp],kind=sp)
        b = cmplx([2.0_sp,4.0_sp,5.0_sp,4.0_sp], &
                  [1.0_sp,1.5_sp,2.0_sp,1.0_sp],kind=sp)
        w = 1.0_sp ! Uniform weights = OLS (weights are always real)
        A_copy = A

        x_weighted = weighted_lstsq(w,A,b,err=state)

        call check(error,state%ok(),'complex weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error,size(x_weighted) == n,'complex weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'complex lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_weighted - x_ols) < tol), &
                   'complex weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_c
    !> Test weighted least-squares with complex data: uniform weights must match lstsq
    subroutine test_weighted_lstsq_z(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: A(m,n),A_copy(m,n),b(m)
        real(dp) :: w(m)
        complex(dp),allocatable :: x_weighted(:),x_ols(:)

        ! Complex test case
        A(:,1) = cmplx(1.0_dp,0.5_dp,kind=dp)
        A(:,2) = cmplx([1.0_dp,2.0_dp,3.0_dp,4.0_dp], &
                       [0.1_dp,0.2_dp,0.3_dp,0.4_dp],kind=dp)
        b = cmplx([2.0_dp,4.0_dp,5.0_dp,4.0_dp], &
                  [1.0_dp,1.5_dp,2.0_dp,1.0_dp],kind=dp)
        w = 1.0_dp ! Uniform weights = OLS (weights are always real)
        A_copy = A

        x_weighted = weighted_lstsq(w,A,b,err=state)

        call check(error,state%ok(),'complex weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error,size(x_weighted) == n,'complex weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'complex lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_weighted - x_ols) < tol), &
                   'complex weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_z
    !> Test weighted least-squares with complex data: uniform weights must match lstsq
    subroutine test_weighted_lstsq_w(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: A(m,n),A_copy(m,n),b(m)
        real(qp) :: w(m)
        complex(qp),allocatable :: x_weighted(:),x_ols(:)

        ! Complex test case
        A(:,1) = cmplx(1.0_qp,0.5_qp,kind=qp)
        A(:,2) = cmplx([1.0_qp,2.0_qp,3.0_qp,4.0_qp], &
                       [0.1_qp,0.2_qp,0.3_qp,0.4_qp],kind=qp)
        b = cmplx([2.0_qp,4.0_qp,5.0_qp,4.0_qp], &
                  [1.0_qp,1.5_qp,2.0_qp,1.0_qp],kind=qp)
        w = 1.0_qp ! Uniform weights = OLS (weights are always real)
        A_copy = A

        x_weighted = weighted_lstsq(w,A,b,err=state)

        call check(error,state%ok(),'complex weighted_lstsq failed: '//state%print())
        if (allocated(error)) return
        
        call check(error,size(x_weighted) == n,'complex weighted_lstsq: wrong solution size')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'complex lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_weighted - x_ols) < tol), &
                   'complex weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_weighted_lstsq_w

    !-------------------------------------------------------------
    !-----     Weighted Least-Squares Subroutine Tests        -----
    !-------------------------------------------------------------

    !> Test basic solve_weighted_lstsq: uniform weights must match ordinary least-squares
    subroutine test_solve_weighted_lstsq_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n),A_original(m,n),A_copy(m,n),b(m),x_sub(n)
        real(sp) :: w(m)
        real(sp),allocatable :: x_ols(:)

        ! Simple test case
        A(:,1) = 1.0_sp
        A(:,2) = [1.0_sp,2.0_sp,3.0_sp,4.0_sp]
        b = [2.0_sp,4.0_sp,5.0_sp,4.0_sp]
        w = 1.0_sp ! Uniform weights = OLS
        A_original = A
        A_copy = A

        call solve_weighted_lstsq(w,A,b,x_sub,err=state)

        call check(error,state%ok(),'solve_weighted_lstsq failed: '//state%print())
        if (allocated(error)) return

        ! Verify A is preserved by default (overwrite_a contract)
        call check(error,all(A == A_original), &
                   'solve_weighted_lstsq must preserve A by default')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares exactly
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_sub - x_ols) < tol), &
                   'solve_weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_s

    !> Test that non-uniform weights change the solution (subroutine version)
    subroutine test_solve_weighted_lstsq_effect_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: A(m,n),b(m),x_uniform(n),x_weighted(n)
        real(sp) :: w_uniform(m),w_nonuniform(m)

        ! Setup problem
        A(:,1) = 1.0_sp
        A(:,2) = [1.0_sp,2.0_sp,3.0_sp,4.0_sp]
        b = [1.0_sp,3.0_sp,2.0_sp,5.0_sp]

        w_uniform = 1.0_sp
        w_nonuniform = [10.0_sp,1.0_sp,1.0_sp,10.0_sp] ! Weight first and last more

        call solve_weighted_lstsq(w_uniform,A,b,x_uniform,err=state)
        call check(error,state%ok(),'uniform solve_weighted_lstsq failed')
        if (allocated(error)) return

        call solve_weighted_lstsq(w_nonuniform,A,b,x_weighted,err=state)
        call check(error,state%ok(),'non-uniform solve_weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error,any(abs(x_uniform - x_weighted) > tol), &
                   'solve_weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_effect_s

    !> Test error on negative weights (subroutine version)
    subroutine test_solve_weighted_lstsq_negative_s(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(sp) :: A(3,2),b(3),x(2)
        real(sp) :: w(3)

        A = 1.0_sp
        b = 1.0_sp
        w = [-1.0_sp,1.0_sp,1.0_sp] ! Invalid: negative weight!

        call solve_weighted_lstsq(w,A,b,x,err=state)

        call check(error,state%error(),'solve_weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_negative_s

    !> Test basic solve_weighted_lstsq: uniform weights must match ordinary least-squares
    subroutine test_solve_weighted_lstsq_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n),A_original(m,n),A_copy(m,n),b(m),x_sub(n)
        real(dp) :: w(m)
        real(dp),allocatable :: x_ols(:)

        ! Simple test case
        A(:,1) = 1.0_dp
        A(:,2) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        b = [2.0_dp,4.0_dp,5.0_dp,4.0_dp]
        w = 1.0_dp ! Uniform weights = OLS
        A_original = A
        A_copy = A

        call solve_weighted_lstsq(w,A,b,x_sub,err=state)

        call check(error,state%ok(),'solve_weighted_lstsq failed: '//state%print())
        if (allocated(error)) return

        ! Verify A is preserved by default (overwrite_a contract)
        call check(error,all(A == A_original), &
                   'solve_weighted_lstsq must preserve A by default')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares exactly
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_sub - x_ols) < tol), &
                   'solve_weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_d

    !> Test that non-uniform weights change the solution (subroutine version)
    subroutine test_solve_weighted_lstsq_effect_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: A(m,n),b(m),x_uniform(n),x_weighted(n)
        real(dp) :: w_uniform(m),w_nonuniform(m)

        ! Setup problem
        A(:,1) = 1.0_dp
        A(:,2) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        b = [1.0_dp,3.0_dp,2.0_dp,5.0_dp]

        w_uniform = 1.0_dp
        w_nonuniform = [10.0_dp,1.0_dp,1.0_dp,10.0_dp] ! Weight first and last more

        call solve_weighted_lstsq(w_uniform,A,b,x_uniform,err=state)
        call check(error,state%ok(),'uniform solve_weighted_lstsq failed')
        if (allocated(error)) return

        call solve_weighted_lstsq(w_nonuniform,A,b,x_weighted,err=state)
        call check(error,state%ok(),'non-uniform solve_weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error,any(abs(x_uniform - x_weighted) > tol), &
                   'solve_weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_effect_d

    !> Test error on negative weights (subroutine version)
    subroutine test_solve_weighted_lstsq_negative_d(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(dp) :: A(3,2),b(3),x(2)
        real(dp) :: w(3)

        A = 1.0_dp
        b = 1.0_dp
        w = [-1.0_dp,1.0_dp,1.0_dp] ! Invalid: negative weight!

        call solve_weighted_lstsq(w,A,b,x,err=state)

        call check(error,state%error(),'solve_weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_negative_d

    !> Test basic solve_weighted_lstsq: uniform weights must match ordinary least-squares
    subroutine test_solve_weighted_lstsq_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: A(m,n),A_original(m,n),A_copy(m,n),b(m),x_sub(n)
        real(qp) :: w(m)
        real(qp),allocatable :: x_ols(:)

        ! Simple test case
        A(:,1) = 1.0_qp
        A(:,2) = [1.0_qp,2.0_qp,3.0_qp,4.0_qp]
        b = [2.0_qp,4.0_qp,5.0_qp,4.0_qp]
        w = 1.0_qp ! Uniform weights = OLS
        A_original = A
        A_copy = A

        call solve_weighted_lstsq(w,A,b,x_sub,err=state)

        call check(error,state%ok(),'solve_weighted_lstsq failed: '//state%print())
        if (allocated(error)) return

        ! Verify A is preserved by default (overwrite_a contract)
        call check(error,all(A == A_original), &
                   'solve_weighted_lstsq must preserve A by default')
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares exactly
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_sub - x_ols) < tol), &
                   'solve_weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_q

    !> Test that non-uniform weights change the solution (subroutine version)
    subroutine test_solve_weighted_lstsq_effect_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: A(m,n),b(m),x_uniform(n),x_weighted(n)
        real(qp) :: w_uniform(m),w_nonuniform(m)

        ! Setup problem
        A(:,1) = 1.0_qp
        A(:,2) = [1.0_qp,2.0_qp,3.0_qp,4.0_qp]
        b = [1.0_qp,3.0_qp,2.0_qp,5.0_qp]

        w_uniform = 1.0_qp
        w_nonuniform = [10.0_qp,1.0_qp,1.0_qp,10.0_qp] ! Weight first and last more

        call solve_weighted_lstsq(w_uniform,A,b,x_uniform,err=state)
        call check(error,state%ok(),'uniform solve_weighted_lstsq failed')
        if (allocated(error)) return

        call solve_weighted_lstsq(w_nonuniform,A,b,x_weighted,err=state)
        call check(error,state%ok(),'non-uniform solve_weighted_lstsq failed')
        if (allocated(error)) return

        ! Solutions should be different
        call check(error,any(abs(x_uniform - x_weighted) > tol), &
                   'solve_weighted_lstsq: non-uniform weights should change solution')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_effect_q

    !> Test error on negative weights (subroutine version)
    subroutine test_solve_weighted_lstsq_negative_q(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        real(qp) :: A(3,2),b(3),x(2)
        real(qp) :: w(3)

        A = 1.0_qp
        b = 1.0_qp
        w = [-1.0_qp,1.0_qp,1.0_qp] ! Invalid: negative weight!

        call solve_weighted_lstsq(w,A,b,x,err=state)

        call check(error,state%error(),'solve_weighted_lstsq should fail on negative weights')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_negative_q

    !> Test solve_weighted_lstsq with complex data: uniform weights must match lstsq
    subroutine test_solve_weighted_lstsq_c(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: A(m,n),A_copy(m,n),b(m),x_sub(n)
        real(sp) :: w(m)
        complex(sp),allocatable :: x_ols(:)

        ! Complex test case
        A(:,1) = cmplx(1.0_sp,0.5_sp,kind=sp)
        A(:,2) = cmplx([1.0_sp,2.0_sp,3.0_sp,4.0_sp], &
                       [0.1_sp,0.2_sp,0.3_sp,0.4_sp],kind=sp)
        b = cmplx([2.0_sp,4.0_sp,5.0_sp,4.0_sp], &
                  [1.0_sp,1.5_sp,2.0_sp,1.0_sp],kind=sp)
        w = 1.0_sp ! Uniform weights = OLS (weights are always real)
        A_copy = A

        call solve_weighted_lstsq(w,A,b,x_sub,err=state)

        call check(error,state%ok(),'complex solve_weighted_lstsq failed: '//state%print())
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'complex lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_sub - x_ols) < tol), &
                   'complex solve_weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_c
    !> Test solve_weighted_lstsq with complex data: uniform weights must match lstsq
    subroutine test_solve_weighted_lstsq_z(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: A(m,n),A_copy(m,n),b(m),x_sub(n)
        real(dp) :: w(m)
        complex(dp),allocatable :: x_ols(:)

        ! Complex test case
        A(:,1) = cmplx(1.0_dp,0.5_dp,kind=dp)
        A(:,2) = cmplx([1.0_dp,2.0_dp,3.0_dp,4.0_dp], &
                       [0.1_dp,0.2_dp,0.3_dp,0.4_dp],kind=dp)
        b = cmplx([2.0_dp,4.0_dp,5.0_dp,4.0_dp], &
                  [1.0_dp,1.5_dp,2.0_dp,1.0_dp],kind=dp)
        w = 1.0_dp ! Uniform weights = OLS (weights are always real)
        A_copy = A

        call solve_weighted_lstsq(w,A,b,x_sub,err=state)

        call check(error,state%ok(),'complex solve_weighted_lstsq failed: '//state%print())
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'complex lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_sub - x_ols) < tol), &
                   'complex solve_weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_z
    !> Test solve_weighted_lstsq with complex data: uniform weights must match lstsq
    subroutine test_solve_weighted_lstsq_w(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state
        integer(ilp),parameter :: m = 4,n = 2
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: A(m,n),A_copy(m,n),b(m),x_sub(n)
        real(qp) :: w(m)
        complex(qp),allocatable :: x_ols(:)

        ! Complex test case
        A(:,1) = cmplx(1.0_qp,0.5_qp,kind=qp)
        A(:,2) = cmplx([1.0_qp,2.0_qp,3.0_qp,4.0_qp], &
                       [0.1_qp,0.2_qp,0.3_qp,0.4_qp],kind=qp)
        b = cmplx([2.0_qp,4.0_qp,5.0_qp,4.0_qp], &
                  [1.0_qp,1.5_qp,2.0_qp,1.0_qp],kind=qp)
        w = 1.0_qp ! Uniform weights = OLS (weights are always real)
        A_copy = A

        call solve_weighted_lstsq(w,A,b,x_sub,err=state)

        call check(error,state%ok(),'complex solve_weighted_lstsq failed: '//state%print())
        if (allocated(error)) return

        ! KEY TEST: Uniform weights should match ordinary least squares
        x_ols = lstsq(A_copy,b,err=state)

        call check(error,state%ok(),'complex lstsq failed: '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(x_sub - x_ols) < tol), &
                   'complex solve_weighted_lstsq with uniform weights must match lstsq')
        if (allocated(error)) return

    end subroutine test_solve_weighted_lstsq_w

    ! Test issue #823
    subroutine test_issue_823(error)
        type(error_type),allocatable,intent(out) :: error
        
        ! Dimension of the problem.
        integer(ilp),parameter :: n = 42
        ! Data for the least-squares problem.
        complex(dp) :: A(n + 1,n),b(n + 1),x_true(n),x_lstsq(n)
        ! Internal variables.
        real(dp),allocatable :: tmp(:,:,:),tmp_vec(:,:)
        ! Error handler
        type(la_state) :: state

        ! Zero-out data.
        A = 0.0_dp
        b = 0.0_dp
        x_lstsq = 0.0_dp
        allocate (tmp(n + 1,n,2),tmp_vec(n,2),source=0.0_dp)

        ! Generate a random complex least-squares problem of size (n+1, n).
        call random_number(tmp)
        call random_number(tmp_vec)
        
        A = cmplx(tmp(:,:,1),tmp(:,:,2),kind=dp)
        x_true = cmplx(tmp_vec(:,1),tmp_vec(:,2),kind=dp)
        b = matmul(A,x_true)

        ! Solve the lstsq problem.
        call solve_lstsq(A,b,x_lstsq,err=state)
          
        ! Check that no segfault occurred
        call check(error,state%ok(),'issue 823 returned '//state%print())
        if (allocated(error)) return

        ! Check that least squares are verified
        call check(error,all(abs(x_true - x_lstsq) < sqrt(epsilon(0.0_dp))),'issue 823 results')
        if (allocated(error)) return

    end subroutine test_issue_823

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

end module test_linalg_least_squares

program test_lstsq
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_least_squares,only:test_least_squares
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_least_squares",test_least_squares) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_lstsq
