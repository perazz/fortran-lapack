! The high-level procedures that are pure, called from a pure procedure: a change that makes
! one of them impure fails to compile this test
module test_linalg_pure
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:eigvals,eigvalsh,eig,svd,invert,schur,solve_lstsq, &
                              generalized_lstsq,solve_generalized_lstsq, &
                              operator(.inv.),operator(.pinv.)
    use la_state_type,only:la_state

    implicit none(type,external)
    private

    public :: test_pure_procedures

    contains

    !> Pure procedure tests
    subroutine test_pure_procedures(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("pure_calls_s",test_pure_calls_s))
        call add_test(tests,new_unittest("pure_calls_d",test_pure_calls_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("pure_calls_x",test_pure_calls_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("pure_calls_q",test_pure_calls_q))
#endif
        call add_test(tests,new_unittest("pure_calls_c",test_pure_calls_c))
        call add_test(tests,new_unittest("pure_calls_z",test_pure_calls_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("pure_calls_y",test_pure_calls_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("pure_calls_w",test_pure_calls_w))
#endif

    end subroutine test_pure_procedures

    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp),parameter :: d(n) = [1.0_sp,2.0_sp,4.0_sp]
        real(sp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(sp) :: lambda(n),lambda_eig(n)
        real(sp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_s(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_s

    pure subroutine pure_calls_s(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        real(sp),intent(in) :: a(:,:)
        complex(sp),intent(out) :: lambda(:),lambda_eig(:)
        real(sp),intent(out) :: lambdah(:),s(:)
        real(sp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        real(sp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_s

    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp),parameter :: d(n) = [1.0_dp,2.0_dp,4.0_dp]
        real(dp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(dp) :: lambda(n),lambda_eig(n)
        real(dp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_d(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_d

    pure subroutine pure_calls_d(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        real(dp),intent(in) :: a(:,:)
        complex(dp),intent(out) :: lambda(:),lambda_eig(:)
        real(dp),intent(out) :: lambdah(:),s(:)
        real(dp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        real(dp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_d

#ifdef LA_WITH_XDP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp),parameter :: d(n) = [1.0_xdp,2.0_xdp,4.0_xdp]
        real(xdp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(xdp) :: lambda(n),lambda_eig(n)
        real(xdp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_x(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_x

    pure subroutine pure_calls_x(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        real(xdp),intent(in) :: a(:,:)
        complex(xdp),intent(out) :: lambda(:),lambda_eig(:)
        real(xdp),intent(out) :: lambdah(:),s(:)
        real(xdp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        real(xdp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_x

#endif
#ifdef LA_WITH_QP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp),parameter :: d(n) = [1.0_qp,2.0_qp,4.0_qp]
        real(qp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(qp) :: lambda(n),lambda_eig(n)
        real(qp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_q(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_q

    pure subroutine pure_calls_q(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        real(qp),intent(in) :: a(:,:)
        complex(qp),intent(out) :: lambda(:),lambda_eig(:)
        real(qp),intent(out) :: lambdah(:),s(:)
        real(qp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        real(qp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_q

#endif
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp),parameter :: d(n) = [1.0_sp,2.0_sp,4.0_sp]
        complex(sp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(sp) :: lambda(n),lambda_eig(n)
        real(sp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_c(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_c

    pure subroutine pure_calls_c(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        complex(sp),intent(in) :: a(:,:)
        complex(sp),intent(out) :: lambda(:),lambda_eig(:)
        real(sp),intent(out) :: lambdah(:),s(:)
        complex(sp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        complex(sp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_c

    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp),parameter :: d(n) = [1.0_dp,2.0_dp,4.0_dp]
        complex(dp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(dp) :: lambda(n),lambda_eig(n)
        real(dp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_z(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_z

    pure subroutine pure_calls_z(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        complex(dp),intent(in) :: a(:,:)
        complex(dp),intent(out) :: lambda(:),lambda_eig(:)
        real(dp),intent(out) :: lambdah(:),s(:)
        complex(dp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        complex(dp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_z

#ifdef LA_WITH_XDP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp),parameter :: d(n) = [1.0_xdp,2.0_xdp,4.0_xdp]
        complex(xdp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(xdp) :: lambda(n),lambda_eig(n)
        real(xdp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_y(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_y

    pure subroutine pure_calls_y(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        complex(xdp),intent(in) :: a(:,:)
        complex(xdp),intent(out) :: lambda(:),lambda_eig(:)
        real(xdp),intent(out) :: lambdah(:),s(:)
        complex(xdp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        complex(xdp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_y

#endif
#ifdef LA_WITH_QP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp),parameter :: d(n) = [1.0_qp,2.0_qp,4.0_qp]
        complex(qp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n)
        complex(qp) :: lambda(n),lambda_eig(n)
        real(qp) :: lambdah(n),s(n)
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_w(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)

        call check(error,ok,'no error state from the pure subroutines')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda - d(i)) < tol),i=1,n)]),'eigvals')
        if (allocated(error)) return

        call check(error,all([(any(abs(lambda_eig - d(i)) < tol),i=1,n)]),'eig')
        if (allocated(error)) return

        call check(error,all(abs(lambdah - d) < tol),'eigvalsh')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svd')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'.inv.')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'.pinv.')
        if (allocated(error)) return

        call check(error,all(abs(inverted - ainv) < tol),'invert')
        if (allocated(error)) return

        call check(error,abs(t(1,1) + t(2,2) + t(3,3) - sum(d)) < tol .and. &
                          abs(t(2,1)) + abs(t(3,1)) + abs(t(3,2)) < tol,'schur')
        if (allocated(error)) return

        call check(error,all(abs(x - 1) < tol),'solve_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xg - 1) < tol),'generalized_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xgs - 1) < tol),'solve_generalized_lstsq')
        if (allocated(error)) return

    end subroutine test_pure_calls_w

    pure subroutine pure_calls_w(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,ok)
        complex(qp),intent(in) :: a(:,:)
        complex(qp),intent(out) :: lambda(:),lambda_eig(:)
        real(qp),intent(out) :: lambdah(:),s(:)
        complex(qp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:)
        logical(lk),intent(out) :: ok

        complex(qp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2))
        type(la_state) :: err(6)
        integer :: i

        lambda = eigvals(a)
        lambdah = eigvalsh(a)
        inva = .inv.a
        pinva = .pinv.a

        work = a
        call eig(work,lambda_eig,err=err(1))

        work = a
        call svd(work,s,err=err(2))

        inverted = a
        call invert(inverted,err=err(3))

        work = a
        call schur(work,t,err=err(4))

        work = a
        call solve_lstsq(work, [(a(i,i),i=1,size(a,1))],x,err=err(5))

        ! a is also a valid covariance: the square system has the exact solution x = 1
        xg = generalized_lstsq(a,a, [(a(i,i),i=1,size(a,1))])

        work = a
        w = a
        call solve_generalized_lstsq(w,work, [(a(i,i),i=1,size(a,1))],xgs,err=err(6))

        ok = all(err%ok())

    end subroutine pure_calls_w

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

end module test_linalg_pure

program test_pure
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_pure,only:test_pure_procedures
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_pure",test_pure_procedures) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_pure
