! The high-level procedures that are pure, called from a pure procedure: a change that makes
! one of them impure fails to compile this test. The function forms with err return the error
! state on a failing input instead of stopping.
module test_linalg_pure
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:eigvals,eigvalsh,eig,svd,invert,schur,solve_lstsq, &
                              generalized_lstsq,solve_generalized_lstsq,qr, &
                              matrix_rank,cond,slogdet,matrix_power, &
                              solve,inv,det,chol,lstsq,weighted_lstsq,constrained_lstsq, &
                              pinv,svdvals,mnorm, &
                              operator(.inv.),operator(.pinv.),operator(.det.)
    use la_state_type,only:la_state,LINALG_ERROR,LINALG_VALUE_ERROR

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
        call add_test(tests,new_unittest("pure_functions_s",test_pure_functions_s))
        call add_test(tests,new_unittest("err_functions_s",test_err_functions_s))
        call add_test(tests,new_unittest("pure_calls_d",test_pure_calls_d))
        call add_test(tests,new_unittest("pure_functions_d",test_pure_functions_d))
        call add_test(tests,new_unittest("err_functions_d",test_err_functions_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("pure_calls_x",test_pure_calls_x))
        call add_test(tests,new_unittest("pure_functions_x",test_pure_functions_x))
        call add_test(tests,new_unittest("err_functions_x",test_err_functions_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("pure_calls_q",test_pure_calls_q))
        call add_test(tests,new_unittest("pure_functions_q",test_pure_functions_q))
        call add_test(tests,new_unittest("err_functions_q",test_err_functions_q))
#endif
        call add_test(tests,new_unittest("pure_calls_c",test_pure_calls_c))
        call add_test(tests,new_unittest("pure_functions_c",test_pure_functions_c))
        call add_test(tests,new_unittest("err_functions_c",test_err_functions_c))
        call add_test(tests,new_unittest("pure_calls_z",test_pure_calls_z))
        call add_test(tests,new_unittest("pure_functions_z",test_pure_functions_z))
        call add_test(tests,new_unittest("err_functions_z",test_err_functions_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("pure_calls_y",test_pure_calls_y))
        call add_test(tests,new_unittest("pure_functions_y",test_pure_functions_y))
        call add_test(tests,new_unittest("err_functions_y",test_err_functions_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("pure_calls_w",test_pure_calls_w))
        call add_test(tests,new_unittest("pure_functions_w",test_pure_functions_w))
        call add_test(tests,new_unittest("err_functions_w",test_err_functions_w))
#endif

    end subroutine test_pure_procedures

    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp),parameter :: d(n) = [1.0_sp,2.0_sp,4.0_sp]
        real(sp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(sp) :: lambda(n),lambda_eig(n)
        real(sp) :: lambdah(n),s(n),c(3),lad
        real(sp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_s(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_s(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_sp,4.0_sp,5.25_sp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_sp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_s

    pure subroutine pure_calls_s(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        real(sp),intent(in) :: a(:,:)
        complex(sp),intent(out) :: lambda(:),lambda_eig(:)
        real(sp),intent(out) :: lambdah(:),s(:)
        real(sp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        real(sp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_s

    pure subroutine pure_scalars_s(a,rnk,c,sgn,lad,a2,am1)
        real(sp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(sp),intent(out) :: c(3),lad
        real(sp),intent(out) :: sgn,a2(:,:),am1(:,:)

        real(sp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_s

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp),parameter :: d(n) = [1.0_sp,2.0_sp,4.0_sp]
        real(sp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(sp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_s(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_sp,sqrt(21.0_sp),sqrt(21.0_sp),sqrt(21.0_sp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_s

    pure subroutine pure_functions_s(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        real(sp),intent(in) :: a(:,:)
        real(sp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(sp),intent(out) :: s(:),nrm(4)

        real(sp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(sp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_s

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        real(sp),allocatable :: x(:),xm(:,:)
        real(sp) :: w(n),nrm
        real(sp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_s

    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp),parameter :: d(n) = [1.0_dp,2.0_dp,4.0_dp]
        real(dp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(dp) :: lambda(n),lambda_eig(n)
        real(dp) :: lambdah(n),s(n),c(3),lad
        real(dp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_d(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_d(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_dp,4.0_dp,5.25_dp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_dp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_d

    pure subroutine pure_calls_d(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        real(dp),intent(in) :: a(:,:)
        complex(dp),intent(out) :: lambda(:),lambda_eig(:)
        real(dp),intent(out) :: lambdah(:),s(:)
        real(dp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        real(dp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_d

    pure subroutine pure_scalars_d(a,rnk,c,sgn,lad,a2,am1)
        real(dp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(dp),intent(out) :: c(3),lad
        real(dp),intent(out) :: sgn,a2(:,:),am1(:,:)

        real(dp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_d

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp),parameter :: d(n) = [1.0_dp,2.0_dp,4.0_dp]
        real(dp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(dp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_d(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_dp,sqrt(21.0_dp),sqrt(21.0_dp),sqrt(21.0_dp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_d

    pure subroutine pure_functions_d(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        real(dp),intent(in) :: a(:,:)
        real(dp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(dp),intent(out) :: s(:),nrm(4)

        real(dp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(dp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_d

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        real(dp),allocatable :: x(:),xm(:,:)
        real(dp) :: w(n),nrm
        real(dp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_d

#ifdef LA_WITH_XDP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp),parameter :: d(n) = [1.0_xdp,2.0_xdp,4.0_xdp]
        real(xdp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(xdp) :: lambda(n),lambda_eig(n)
        real(xdp) :: lambdah(n),s(n),c(3),lad
        real(xdp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_x(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_x(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_xdp,4.0_xdp,5.25_xdp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_xdp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_x

    pure subroutine pure_calls_x(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        real(xdp),intent(in) :: a(:,:)
        complex(xdp),intent(out) :: lambda(:),lambda_eig(:)
        real(xdp),intent(out) :: lambdah(:),s(:)
        real(xdp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        real(xdp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_x

    pure subroutine pure_scalars_x(a,rnk,c,sgn,lad,a2,am1)
        real(xdp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(xdp),intent(out) :: c(3),lad
        real(xdp),intent(out) :: sgn,a2(:,:),am1(:,:)

        real(xdp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_x

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp),parameter :: d(n) = [1.0_xdp,2.0_xdp,4.0_xdp]
        real(xdp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(xdp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_x(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_xdp,sqrt(21.0_xdp),sqrt(21.0_xdp),sqrt(21.0_xdp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_x

    pure subroutine pure_functions_x(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        real(xdp),intent(in) :: a(:,:)
        real(xdp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(xdp),intent(out) :: s(:),nrm(4)

        real(xdp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(xdp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_x

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        real(xdp),allocatable :: x(:),xm(:,:)
        real(xdp) :: w(n),nrm
        real(xdp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_x

#endif
#ifdef LA_WITH_QP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp),parameter :: d(n) = [1.0_qp,2.0_qp,4.0_qp]
        real(qp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(qp) :: lambda(n),lambda_eig(n)
        real(qp) :: lambdah(n),s(n),c(3),lad
        real(qp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_q(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_q(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_qp,4.0_qp,5.25_qp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_qp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_q

    pure subroutine pure_calls_q(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        real(qp),intent(in) :: a(:,:)
        complex(qp),intent(out) :: lambda(:),lambda_eig(:)
        real(qp),intent(out) :: lambdah(:),s(:)
        real(qp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        real(qp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_q

    pure subroutine pure_scalars_q(a,rnk,c,sgn,lad,a2,am1)
        real(qp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(qp),intent(out) :: c(3),lad
        real(qp),intent(out) :: sgn,a2(:,:),am1(:,:)

        real(qp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_q

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp),parameter :: d(n) = [1.0_qp,2.0_qp,4.0_qp]
        real(qp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(qp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_q(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_qp,sqrt(21.0_qp),sqrt(21.0_qp),sqrt(21.0_qp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_q

    pure subroutine pure_functions_q(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        real(qp),intent(in) :: a(:,:)
        real(qp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(qp),intent(out) :: s(:),nrm(4)

        real(qp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(qp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_q

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        real(qp),allocatable :: x(:),xm(:,:)
        real(qp) :: w(n),nrm
        real(qp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_q

#endif
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp),parameter :: d(n) = [1.0_sp,2.0_sp,4.0_sp]
        complex(sp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(sp) :: lambda(n),lambda_eig(n)
        real(sp) :: lambdah(n),s(n),c(3),lad
        complex(sp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_c(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_c(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_sp,4.0_sp,5.25_sp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_sp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_c

    pure subroutine pure_calls_c(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        complex(sp),intent(in) :: a(:,:)
        complex(sp),intent(out) :: lambda(:),lambda_eig(:)
        real(sp),intent(out) :: lambdah(:),s(:)
        complex(sp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        complex(sp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_c

    pure subroutine pure_scalars_c(a,rnk,c,sgn,lad,a2,am1)
        complex(sp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(sp),intent(out) :: c(3),lad
        complex(sp),intent(out) :: sgn,a2(:,:),am1(:,:)

        complex(sp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_c

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp),parameter :: d(n) = [1.0_sp,2.0_sp,4.0_sp]
        complex(sp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(sp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_c(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_sp,sqrt(21.0_sp),sqrt(21.0_sp),sqrt(21.0_sp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_c

    pure subroutine pure_functions_c(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        complex(sp),intent(in) :: a(:,:)
        complex(sp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(sp),intent(out) :: s(:),nrm(4)

        complex(sp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(sp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_c

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(sp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        complex(sp),allocatable :: x(:),xm(:,:)
        real(sp) :: w(n),nrm
        real(sp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_c

    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp),parameter :: d(n) = [1.0_dp,2.0_dp,4.0_dp]
        complex(dp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(dp) :: lambda(n),lambda_eig(n)
        real(dp) :: lambdah(n),s(n),c(3),lad
        complex(dp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_z(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_z(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_dp,4.0_dp,5.25_dp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_dp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_z

    pure subroutine pure_calls_z(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        complex(dp),intent(in) :: a(:,:)
        complex(dp),intent(out) :: lambda(:),lambda_eig(:)
        real(dp),intent(out) :: lambdah(:),s(:)
        complex(dp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        complex(dp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_z

    pure subroutine pure_scalars_z(a,rnk,c,sgn,lad,a2,am1)
        complex(dp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(dp),intent(out) :: c(3),lad
        complex(dp),intent(out) :: sgn,a2(:,:),am1(:,:)

        complex(dp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_z

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp),parameter :: d(n) = [1.0_dp,2.0_dp,4.0_dp]
        complex(dp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(dp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_z(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_dp,sqrt(21.0_dp),sqrt(21.0_dp),sqrt(21.0_dp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_z

    pure subroutine pure_functions_z(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        complex(dp),intent(in) :: a(:,:)
        complex(dp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(dp),intent(out) :: s(:),nrm(4)

        complex(dp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(dp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_z

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(dp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        complex(dp),allocatable :: x(:),xm(:,:)
        real(dp) :: w(n),nrm
        real(dp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_z

#ifdef LA_WITH_XDP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp),parameter :: d(n) = [1.0_xdp,2.0_xdp,4.0_xdp]
        complex(xdp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(xdp) :: lambda(n),lambda_eig(n)
        real(xdp) :: lambdah(n),s(n),c(3),lad
        complex(xdp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_y(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_y(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_xdp,4.0_xdp,5.25_xdp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_xdp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_y

    pure subroutine pure_calls_y(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        complex(xdp),intent(in) :: a(:,:)
        complex(xdp),intent(out) :: lambda(:),lambda_eig(:)
        real(xdp),intent(out) :: lambdah(:),s(:)
        complex(xdp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        complex(xdp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_y

    pure subroutine pure_scalars_y(a,rnk,c,sgn,lad,a2,am1)
        complex(xdp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(xdp),intent(out) :: c(3),lad
        complex(xdp),intent(out) :: sgn,a2(:,:),am1(:,:)

        complex(xdp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_y

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp),parameter :: d(n) = [1.0_xdp,2.0_xdp,4.0_xdp]
        complex(xdp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(xdp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_y(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_xdp,sqrt(21.0_xdp),sqrt(21.0_xdp),sqrt(21.0_xdp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_y

    pure subroutine pure_functions_y(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        complex(xdp),intent(in) :: a(:,:)
        complex(xdp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(xdp),intent(out) :: s(:),nrm(4)

        complex(xdp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(xdp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_y

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(xdp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        complex(xdp),allocatable :: x(:),xm(:,:)
        real(xdp) :: w(n),nrm
        real(xdp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_y

#endif
#ifdef LA_WITH_QP
    !> Run every pure procedure on A = diag(1,2,4), whose results are known
    subroutine test_pure_calls_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp),parameter :: d(n) = [1.0_qp,2.0_qp,4.0_qp]
        complex(qp) :: a(n,n),inva(n,n),pinva(n,n),inverted(n,n),t(n,n),x(n),ainv(n,n), &
                  xg(n),xgs(n),r(n,n)
        complex(qp) :: lambda(n),lambda_eig(n)
        real(qp) :: lambdah(n),s(n),c(3),lad
        complex(qp) :: sgn,a2(n,n),am1(n,n)
        integer(ilp) :: pivots(n),rnk
        logical(lk) :: ok
        integer :: i

        a = 0
        ainv = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
        end do

        call pure_calls_w(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        call pure_scalars_w(a,rnk,c,sgn,lad,a2,am1)

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

        call check(error,all(pivots == [3,2,1]) .and. &
                          all(abs([(abs(r(i,i)),i=1,n)] - d(n:1:-1)) < tol),'pivoted qr')
        if (allocated(error)) return

        call check(error,rnk == n,'matrix_rank')
        if (allocated(error)) return

        ! 2- and 1-norms: 4/1; Frobenius norm: sqrt(21)*sqrt(21/16)
        call check(error,all(abs(c - [4.0_qp,4.0_qp,5.25_qp]) < tol),'cond')
        if (allocated(error)) return

        call check(error,abs(sgn - 1) < tol .and. abs(lad - log(8.0_qp)) < tol,'slogdet')
        if (allocated(error)) return

        call check(error,all(abs(a2 - matmul(a,a)) < tol) .and. all(abs(am1 - ainv) < tol),'matrix_power')
        if (allocated(error)) return

    end subroutine test_pure_calls_w

    pure subroutine pure_calls_w(a,lambda,lambda_eig,lambdah,s,inva,pinva,inverted,t,x,xg,xgs,r,pivots,ok)
        complex(qp),intent(in) :: a(:,:)
        complex(qp),intent(out) :: lambda(:),lambda_eig(:)
        real(qp),intent(out) :: lambdah(:),s(:)
        complex(qp),intent(out) :: inva(:,:),pinva(:,:),inverted(:,:),t(:,:),x(:),xg(:),xgs(:),r(:,:)
        integer(ilp),intent(out) :: pivots(:)
        logical(lk),intent(out) :: ok

        complex(qp) :: work(size(a,1),size(a,2)),w(size(a,1),size(a,2)),q(size(a,1),size(a,1))
        type(la_state) :: err(7)
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

        work = a
        call qr(work,q,r,pivots,err=err(7))

        ok = all(err%ok())

    end subroutine pure_calls_w

    pure subroutine pure_scalars_w(a,rnk,c,sgn,lad,a2,am1)
        complex(qp),intent(in) :: a(:,:)
        integer(ilp),intent(out) :: rnk
        real(qp),intent(out) :: c(3),lad
        complex(qp),intent(out) :: sgn,a2(:,:),am1(:,:)

        complex(qp) :: work(size(a,1),size(a,2))

        rnk = matrix_rank(a)
        c = [cond(a),cond(a,1_ilp),cond(a,'fro')]
        a2 = matrix_power(a,2_ilp)
        am1 = matrix_power(a,-1_ilp)

        work = a
        call slogdet(work,sgn,lad)

    end subroutine pure_scalars_w

    !> Run every pure function form on A = diag(1,2,4), whose results are known
    subroutine test_pure_functions_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp),parameter :: d(n) = [1.0_qp,2.0_qp,4.0_qp]
        complex(qp) :: a(n,n),ainv(n,n),lref(n,n),x(n),inva(n,n),dets(2),l(n,n),xl(n),xw(n), &
                  xc(n),pinva(n,n)
        real(qp) :: s(n),nrm(4)
        integer :: i

        a = 0
        ainv = 0
        lref = 0
        do i = 1,n
            a(i,i) = d(i)
            ainv(i,i) = 1/d(i)
            lref(i,i) = sqrt(d(i))
        end do

        call pure_functions_w(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)

        call check(error,all(abs(x - 1) < tol),'solve')
        if (allocated(error)) return

        call check(error,all(abs(inva - ainv) < tol),'inv')
        if (allocated(error)) return

        call check(error,all(abs(dets - product(d)) < tol),'det, .det.')
        if (allocated(error)) return

        call check(error,all(abs(l - lref) < tol),'chol')
        if (allocated(error)) return

        call check(error,all(abs(xl - 1) < tol),'lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xw - 1) < tol),'weighted_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(xc - 1) < tol),'constrained_lstsq')
        if (allocated(error)) return

        call check(error,all(abs(pinva - ainv) < tol),'pinv')
        if (allocated(error)) return

        call check(error,all(abs(s - d(n:1:-1)) < tol),'svdvals')
        if (allocated(error)) return

        ! 1-norm: 4; Frobenius norm: sqrt(21), of the matrix and of both slices of a rank-3 array
        call check(error,all(abs(nrm - [4.0_qp,sqrt(21.0_qp),sqrt(21.0_qp),sqrt(21.0_qp)]) < tol), &
                   'mnorm')
        if (allocated(error)) return

    end subroutine test_pure_functions_w

    pure subroutine pure_functions_w(a,x,inva,dets,l,xl,xw,xc,pinva,s,nrm)
        complex(qp),intent(in) :: a(:,:)
        complex(qp),intent(out) :: x(:),inva(:,:),dets(2),l(:,:),xl(:),xw(:),xc(:),pinva(:,:)
        real(qp),intent(out) :: s(:),nrm(4)

        complex(qp) :: b(size(a,1)),c(1,size(a,2)),a3(size(a,1),size(a,2),2)
        real(qp) :: w(size(a,1))
        integer :: i

        b = [(a(i,i),i=1,size(a,1))]
        w = 1
        c = 0
        c(1,1) = 1
        a3(:,:,1) = a
        a3(:,:,2) = a

        x = solve(a,b)
        inva = inv(a)
        dets = [det(a),.det.a]
        l = chol(a)
        xl = lstsq(a,b)
        xw = weighted_lstsq(w,a,b)
        ! The constraint x(1) = 1 holds at the unconstrained solution x = 1
        xc = constrained_lstsq(a,b,c, [b(1)])
        pinva = pinv(a)
        s = svdvals(a)
        nrm(1:2) = [mnorm(a,'1'),mnorm(a,'fro')]
        nrm(3:4) = mnorm(a3,'fro')

    end subroutine pure_functions_w

    !> The function forms with err return the error state, and do not stop, on a failing input
    subroutine test_err_functions_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(qp) :: z(n,n),rect(n,n - 1),b(n),c(1,n),d(1),empty(0,0),a3(n,n,2),dt
        complex(qp),allocatable :: x(:),xm(:,:)
        real(qp) :: w(n),nrm
        real(qp),allocatable :: s(:),nrm3(:)
        type(la_state) :: err

        z = 0
        rect = 1
        b = 1
        c = 1
        d = 1
        w = -1
        a3 = 1

        x = solve(z,b,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(x) == 0,'solve: singular matrix')
        if (allocated(error)) return

        xm = inv(z,err=err)
        call check(error,err%state == LINALG_ERROR .and. size(xm) == 0,'inv: singular matrix')
        if (allocated(error)) return

        dt = det(rect,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'det: non-square matrix')
        if (allocated(error)) return

        xm = chol(z,err=err)
        call check(error,err%state == LINALG_ERROR,'chol: not positive definite')
        if (allocated(error)) return

        x = lstsq(rect,b(1:n - 1),err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(x) == 0,'lstsq: b has too few rows')
        if (allocated(error)) return

        x = weighted_lstsq(w,rect,b,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'weighted_lstsq: negative weights')
        if (allocated(error)) return

        x = constrained_lstsq(rect,b,c,d,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'constrained_lstsq: a and c columns differ')
        if (allocated(error)) return

        xm = pinv(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'pinv: empty matrix')
        if (allocated(error)) return

        s = svdvals(empty,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'svdvals: empty matrix')
        if (allocated(error)) return

        nrm = mnorm(rect,3_ilp,err=err)
        call check(error,err%state == LINALG_VALUE_ERROR,'mnorm: order 3 is not a matrix norm')
        if (allocated(error)) return

        nrm3 = mnorm(a3,'fro',dim=[1_ilp,1_ilp],err=err)
        call check(error,err%state == LINALG_VALUE_ERROR .and. size(nrm3) == 0,'mnorm: repeated dim')
        if (allocated(error)) return

    end subroutine test_err_functions_w

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
