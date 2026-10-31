! Test QR factorization with column pivoting
module test_linalg_pivoting_qr
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use la_state_type,only:LINALG_VALUE_ERROR,la_state
    use linear_algebra,only:qr,qr_space
    use ieee_arithmetic,only:ieee_value,ieee_quiet_nan

    implicit none(type,external)

    public :: test_pivoting_qr_factorization

    contains

    !> QR factorization tests
    subroutine test_pivoting_qr_factorization(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_s",test_pivoting_qr_random_tall_matrix_s))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_s",test_pivoting_qr_random_rank_deficient_s))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_s",test_pivoting_qr_random_wide_matrix_s))
        call add_test(tests,new_unittest("pivoting_qr_invalid_s",test_pivoting_qr_invalid_s))
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_d",test_pivoting_qr_random_tall_matrix_d))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_d",test_pivoting_qr_random_rank_deficient_d))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_d",test_pivoting_qr_random_wide_matrix_d))
        call add_test(tests,new_unittest("pivoting_qr_invalid_d",test_pivoting_qr_invalid_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_x",test_pivoting_qr_random_tall_matrix_x))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_x",test_pivoting_qr_random_rank_deficient_x))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_x",test_pivoting_qr_random_wide_matrix_x))
        call add_test(tests,new_unittest("pivoting_qr_invalid_x",test_pivoting_qr_invalid_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_q",test_pivoting_qr_random_tall_matrix_q))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_q",test_pivoting_qr_random_rank_deficient_q))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_q",test_pivoting_qr_random_wide_matrix_q))
        call add_test(tests,new_unittest("pivoting_qr_invalid_q",test_pivoting_qr_invalid_q))
#endif
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_c",test_pivoting_qr_random_tall_matrix_c))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_c",test_pivoting_qr_random_rank_deficient_c))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_c",test_pivoting_qr_random_wide_matrix_c))
        call add_test(tests,new_unittest("pivoting_qr_invalid_c",test_pivoting_qr_invalid_c))
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_z",test_pivoting_qr_random_tall_matrix_z))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_z",test_pivoting_qr_random_rank_deficient_z))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_z",test_pivoting_qr_random_wide_matrix_z))
        call add_test(tests,new_unittest("pivoting_qr_invalid_z",test_pivoting_qr_invalid_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_y",test_pivoting_qr_random_tall_matrix_y))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_y",test_pivoting_qr_random_rank_deficient_y))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_y",test_pivoting_qr_random_wide_matrix_y))
        call add_test(tests,new_unittest("pivoting_qr_invalid_y",test_pivoting_qr_invalid_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("pivoting_qr_random_tall_matrix_w",test_pivoting_qr_random_tall_matrix_w))
        call add_test(tests,new_unittest("pivoting_qr_random_rank_deficient_w",test_pivoting_qr_random_rank_deficient_w))
        call add_test(tests,new_unittest("pivoting_qr_random_wide_matrix_w",test_pivoting_qr_random_wide_matrix_w))
        call add_test(tests,new_unittest("pivoting_qr_invalid_w",test_pivoting_qr_invalid_w))
#endif

    end subroutine test_pivoting_qr_factorization

    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(sp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(sp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_s(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_s

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(sp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(sp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        a(:,3) = 0.0_sp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_s

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(sp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(sp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_s(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_s

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        real(sp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(sp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_s

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_s(r) result(ok)
        real(sp),intent(in) :: r(:,:)

        real(sp),parameter :: slack = 1 + sqrt(epsilon(0.0_sp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_s

    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(dp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(dp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_d(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_d

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(dp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(dp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        a(:,3) = 0.0_dp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_d

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        real(dp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(dp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(dp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_d(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_d

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        real(dp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(dp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_d

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_d(r) result(ok)
        real(dp),intent(in) :: r(:,:)

        real(dp),parameter :: slack = 1 + sqrt(epsilon(0.0_dp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_d

#ifdef LA_WITH_XDP
    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(xdp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(xdp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_xdp,ieee_quiet_nan)
        r = ieee_value(0.0_xdp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_x(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_x

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(xdp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(xdp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        a(:,3) = 0.0_xdp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_xdp,ieee_quiet_nan)
        r = ieee_value(0.0_xdp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_x

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(xdp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(xdp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_xdp,ieee_quiet_nan)
        r = ieee_value(0.0_xdp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_x(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_x

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        real(xdp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(xdp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_x

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_x(r) result(ok)
        real(xdp),intent(in) :: r(:,:)

        real(xdp),parameter :: slack = 1 + sqrt(epsilon(0.0_xdp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_x
#endif

#ifdef LA_WITH_QP
    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(qp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(qp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_qp,ieee_quiet_nan)
        r = ieee_value(0.0_qp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_q(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_q

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(qp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(qp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        a(:,3) = 0.0_qp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_qp,ieee_quiet_nan)
        r = ieee_value(0.0_qp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_q

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(qp) :: rea(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        real(qp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        a = rea
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_qp,ieee_quiet_nan)
        r = ieee_value(0.0_qp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_q(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_q

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        real(qp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(qp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_q

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_q(r) result(ok)
        real(qp),intent(in) :: r(:,:)

        real(qp),parameter :: slack = 1 + sqrt(epsilon(0.0_qp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_q
#endif

    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(sp) :: rea(m,n)
        real(sp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(sp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_c(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_c

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(sp) :: rea(m,n)
        real(sp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(sp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        a(:,3) = 0.0_sp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_c

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(sp) :: rea(m,n)
        real(sp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(sp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_c(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_c

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        complex(sp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(sp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_c

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_c(r) result(ok)
        complex(sp),intent(in) :: r(:,:)

        real(sp),parameter :: slack = 1 + sqrt(epsilon(0.0_sp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_c

    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(dp) :: rea(m,n)
        real(dp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(dp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_z(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_z

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(dp) :: rea(m,n)
        real(dp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(dp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        a(:,3) = 0.0_dp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_z

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(dp) :: rea(m,n)
        real(dp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(dp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_z(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_z

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        complex(dp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(dp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_z

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_z(r) result(ok)
        complex(dp),intent(in) :: r(:,:)

        real(dp),parameter :: slack = 1 + sqrt(epsilon(0.0_dp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_z

#ifdef LA_WITH_XDP
    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(xdp) :: rea(m,n)
        real(xdp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(xdp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=xdp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_xdp,ieee_quiet_nan)
        r = ieee_value(0.0_xdp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_y(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_y

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(xdp) :: rea(m,n)
        real(xdp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(xdp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=xdp)
        a(:,3) = 0.0_xdp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_xdp,ieee_quiet_nan)
        r = ieee_value(0.0_xdp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_y

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(xdp),parameter :: tol = 100*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(xdp) :: rea(m,n)
        real(xdp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(xdp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=xdp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_xdp,ieee_quiet_nan)
        r = ieee_value(0.0_xdp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_y(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_y

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        complex(xdp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(xdp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_y

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_y(r) result(ok)
        complex(xdp),intent(in) :: r(:,:)

        real(xdp),parameter :: slack = 1 + sqrt(epsilon(0.0_xdp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_y
#endif

#ifdef LA_WITH_QP
    !> Pivoted QR factorization of a random tall matrix
    subroutine test_pivoting_qr_random_tall_matrix_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(qp) :: rea(m,n)
        real(qp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(qp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_qp,ieee_quiet_nan)
        r = ieee_value(0.0_qp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_w(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_tall_matrix_w

    !> Pivoted QR factorization of a random tall matrix with a zero column
    subroutine test_pivoting_qr_random_rank_deficient_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(qp) :: rea(m,n)
        real(qp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(qp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        a(:,3) = 0.0_qp ! Zero-out column to test rank-deficiency.
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_qp,ieee_quiet_nan)
        r = ieee_value(0.0_qp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! The zero column is pivoted last, and leaves a zero at the end of diag(R)
        call check(error,pivots(n) == 3 .and. abs(r(n,n)) < tol,'zero column pivoted last')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_rank_deficient_w

    !> Pivoted QR factorization of a random wide matrix
    subroutine test_pivoting_qr_random_wide_matrix_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        integer(ilp),parameter :: n = 15_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),aorig(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n),qerr(m,k - 1),rerr(k - 1,n)
        real(qp) :: rea(m,n)
        real(qp) :: ima(m,n)
        integer(ilp) :: pivots(n)
        integer(ilp) :: lwork
        complex(qp),allocatable :: work(:)
        type(la_state) :: state

        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        aorig = a

        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = ieee_value(0.0_qp,ieee_quiet_nan)
        r = ieee_value(0.0_qp,ieee_quiet_nan)
        call qr(a,q,r,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (full)')
        if (allocated(error)) return

        ! Check the pivoting: |diag(R)| does not increase
        call check(error,is_non_increasing_w(r),'non-increasing |diag(R)| (full)')
        if (allocated(error)) return

        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,pivots,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(qred,rred)) < tol),'converged solution (reduced)')
        if (allocated(error)) return

        ! 3) overwrite A
        call qr(a,qred,rred,pivots,overwrite_a=.true.,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(aorig(:,pivots) - matmul(qred,rred)) < tol),'converged solution (overwrite A)')
        if (allocated(error)) return

        ! 4) External storage option
        a = aorig
        call qr_space(a,lwork,pivoting=.true.)
        allocate (work(lwork))
        call qr(a,q,r,pivots,storage=work,err=state)

        ! Check return code
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Check solution
        call check(error,all(abs(a(:,pivots) - matmul(q,r)) < tol),'converged solution (external storage)')
        if (allocated(error)) return

        ! Check that an invalid problem size returns an error
        a = aorig
        call qr(a,qerr,rerr,pivots,err=state)
        call check(error,state%error(),'invalid matrix sizes')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_random_wide_matrix_w

    !> Invalid arguments: pivots, storage, q or r of the wrong size
    subroutine test_pivoting_qr_invalid_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 5_ilp
        complex(qp) :: a(m,n),q(m,m),r(m,n),work(1),qtall(m + 1,m),rwide(m,n + 1)
        real(qp) :: rea(m,n)
        integer(ilp) :: pivots(n),pivots_short(n - 1),pivots_long(n + 1),lwork,lwork_qr
        type(la_state) :: state

        call random_number(rea)
        a = rea

        call qr(a,q,r,pivots_short,err=state)
        call check(error,state%error(),'pivots shorter than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots_long,err=state)
        call check(error,state%error(),'pivots longer than n')
        if (allocated(error)) return

        call qr(a,q,r,pivots,storage=work,err=state)
        call check(error,state%error(),'insufficient storage')
        if (allocated(error)) return

        ! Q and R larger than the factorization can fill
        call qr(a,qtall,r,pivots,err=state)
        call check(error,state%error(),'q with more rows than a')
        if (allocated(error)) return

        call qr(a,q,rwide,pivots,err=state)
        call check(error,state%error(),'r with more columns than a')
        if (allocated(error)) return

        call qr(a,qtall,r,err=state)
        call check(error,state%error(),'q with more rows than a, no pivoting')
        if (allocated(error)) return

        ! pivoting=.false. sizes the workspace of the unpivoted factorization
        call qr_space(a,lwork,pivoting=.false.,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call qr_space(a,lwork_qr,err=state)
        call check(error,state%ok() .and. lwork == lwork_qr,'pivoting=.false. workspace')
        if (allocated(error)) return

    end subroutine test_pivoting_qr_invalid_w

    !> |r(i,i)| does not increase along the diagonal, up to roundoff
    pure logical function is_non_increasing_w(r) result(ok)
        complex(qp),intent(in) :: r(:,:)

        real(qp),parameter :: slack = 1 + sqrt(epsilon(0.0_qp))
        integer :: i

        ok = all([(abs(r(i + 1,i + 1)) <= slack*abs(r(i,i)),i=1,min(size(r,1),size(r,2)) - 1)])

    end function is_non_increasing_w
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

end module test_linalg_pivoting_qr

program test_pivoting_qr
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_pivoting_qr,only:test_pivoting_qr_factorization
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_pivoting_qr",test_pivoting_qr_factorization) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_pivoting_qr
