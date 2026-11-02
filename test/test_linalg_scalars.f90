! Test matrix rank, condition number, log-determinant and matrix power
module test_linalg_scalars
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use la_state_type,only:LINALG_VALUE_ERROR,la_state
    use linear_algebra,only:matrix_rank,cond,slogdet,matrix_power,det,inv,mnorm,svdvals
    use ieee_arithmetic,only:ieee_is_finite,ieee_is_negative

    implicit none(type,external)

    public :: test_matrix_scalars

    contains

    !> Matrix rank, condition number, log-determinant and matrix power tests
    subroutine test_matrix_scalars(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("matrix_rank_ladder_s",test_matrix_rank_ladder_s))
        call add_test(tests,new_unittest("matrix_rank_tolerances_s",test_matrix_rank_tolerances_s))
        call add_test(tests,new_unittest("cond_s",test_cond_s))
        call add_test(tests,new_unittest("cond_invalid_s",test_cond_invalid_s))
        call add_test(tests,new_unittest("slogdet_s",test_slogdet_s))
        call add_test(tests,new_unittest("slogdet_special_s",test_slogdet_special_s))
        call add_test(tests,new_unittest("matrix_power_s",test_matrix_power_s))
        call add_test(tests,new_unittest("matrix_power_invalid_s",test_matrix_power_invalid_s))
        call add_test(tests,new_unittest("matrix_rank_ladder_d",test_matrix_rank_ladder_d))
        call add_test(tests,new_unittest("matrix_rank_tolerances_d",test_matrix_rank_tolerances_d))
        call add_test(tests,new_unittest("cond_d",test_cond_d))
        call add_test(tests,new_unittest("cond_invalid_d",test_cond_invalid_d))
        call add_test(tests,new_unittest("slogdet_d",test_slogdet_d))
        call add_test(tests,new_unittest("slogdet_special_d",test_slogdet_special_d))
        call add_test(tests,new_unittest("matrix_power_d",test_matrix_power_d))
        call add_test(tests,new_unittest("matrix_power_invalid_d",test_matrix_power_invalid_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("matrix_rank_ladder_x",test_matrix_rank_ladder_x))
        call add_test(tests,new_unittest("matrix_rank_tolerances_x",test_matrix_rank_tolerances_x))
        call add_test(tests,new_unittest("cond_x",test_cond_x))
        call add_test(tests,new_unittest("cond_invalid_x",test_cond_invalid_x))
        call add_test(tests,new_unittest("slogdet_x",test_slogdet_x))
        call add_test(tests,new_unittest("slogdet_special_x",test_slogdet_special_x))
        call add_test(tests,new_unittest("matrix_power_x",test_matrix_power_x))
        call add_test(tests,new_unittest("matrix_power_invalid_x",test_matrix_power_invalid_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("matrix_rank_ladder_q",test_matrix_rank_ladder_q))
        call add_test(tests,new_unittest("matrix_rank_tolerances_q",test_matrix_rank_tolerances_q))
        call add_test(tests,new_unittest("cond_q",test_cond_q))
        call add_test(tests,new_unittest("cond_invalid_q",test_cond_invalid_q))
        call add_test(tests,new_unittest("slogdet_q",test_slogdet_q))
        call add_test(tests,new_unittest("slogdet_special_q",test_slogdet_special_q))
        call add_test(tests,new_unittest("matrix_power_q",test_matrix_power_q))
        call add_test(tests,new_unittest("matrix_power_invalid_q",test_matrix_power_invalid_q))
#endif
        call add_test(tests,new_unittest("matrix_rank_ladder_c",test_matrix_rank_ladder_c))
        call add_test(tests,new_unittest("matrix_rank_tolerances_c",test_matrix_rank_tolerances_c))
        call add_test(tests,new_unittest("cond_c",test_cond_c))
        call add_test(tests,new_unittest("cond_invalid_c",test_cond_invalid_c))
        call add_test(tests,new_unittest("slogdet_c",test_slogdet_c))
        call add_test(tests,new_unittest("slogdet_special_c",test_slogdet_special_c))
        call add_test(tests,new_unittest("matrix_power_c",test_matrix_power_c))
        call add_test(tests,new_unittest("matrix_power_invalid_c",test_matrix_power_invalid_c))
        call add_test(tests,new_unittest("matrix_rank_ladder_z",test_matrix_rank_ladder_z))
        call add_test(tests,new_unittest("matrix_rank_tolerances_z",test_matrix_rank_tolerances_z))
        call add_test(tests,new_unittest("cond_z",test_cond_z))
        call add_test(tests,new_unittest("cond_invalid_z",test_cond_invalid_z))
        call add_test(tests,new_unittest("slogdet_z",test_slogdet_z))
        call add_test(tests,new_unittest("slogdet_special_z",test_slogdet_special_z))
        call add_test(tests,new_unittest("matrix_power_z",test_matrix_power_z))
        call add_test(tests,new_unittest("matrix_power_invalid_z",test_matrix_power_invalid_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("matrix_rank_ladder_y",test_matrix_rank_ladder_y))
        call add_test(tests,new_unittest("matrix_rank_tolerances_y",test_matrix_rank_tolerances_y))
        call add_test(tests,new_unittest("cond_y",test_cond_y))
        call add_test(tests,new_unittest("cond_invalid_y",test_cond_invalid_y))
        call add_test(tests,new_unittest("slogdet_y",test_slogdet_y))
        call add_test(tests,new_unittest("slogdet_special_y",test_slogdet_special_y))
        call add_test(tests,new_unittest("matrix_power_y",test_matrix_power_y))
        call add_test(tests,new_unittest("matrix_power_invalid_y",test_matrix_power_invalid_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("matrix_rank_ladder_w",test_matrix_rank_ladder_w))
        call add_test(tests,new_unittest("matrix_rank_tolerances_w",test_matrix_rank_tolerances_w))
        call add_test(tests,new_unittest("cond_w",test_cond_w))
        call add_test(tests,new_unittest("cond_invalid_w",test_cond_invalid_w))
        call add_test(tests,new_unittest("slogdet_w",test_slogdet_w))
        call add_test(tests,new_unittest("slogdet_special_w",test_slogdet_special_w))
        call add_test(tests,new_unittest("matrix_power_w",test_matrix_power_w))
        call add_test(tests,new_unittest("matrix_power_invalid_w",test_matrix_power_invalid_w))
#endif

    end subroutine test_matrix_scalars

    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        real(sp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_sp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = i + j
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_s

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_sp
        d(1,1) = 1.0_sp
        d(2,2) = 1.0e-2_sp
        d(3,3) = 1.0e-4_sp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_sp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_sp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_sp,atol=1.0e-1_sp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_sp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_sp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_s

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(sp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_s(n,n,shift=1.0_sp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_sp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_sp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_sp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(sp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(sp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              real(sp) :: d(3,3)
              d = 0.0_sp
              d(1,1) = 1.0_sp
              d(2,2) = -2.0_sp
              d(3,3) = 4.0_sp
           end function diagonal

    end subroutine test_cond_s

    !> Inputs that cond rejects
    subroutine test_cond_invalid_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: a(3,3),rect(3,2),empty(0,0)
        real(sp) :: c
        type(la_state) :: state

        a = 1.0_sp
        rect = 1.0_sp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_s

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp),parameter :: tol_hilbert = 1.0e5_sp*epsilon(0.0_sp)
        real(sp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(sp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_s(n,n,shift=2.0_sp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_sp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(sp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=sp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_s

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(sp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative entry
        d3 = 0.0_sp
        d3(1,1) = 1.0_sp
        d3(2,2) = -2.0_sp
        d3(3,3) = 3.0_sp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad - log(6.0_sp)) <= tol,'negative determinant')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_sp
        do i = 1,50
           big(i,i) = 10.0_sp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_sp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_sp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_sp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_s

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(sp),parameter :: tol = 1000*epsilon(0.0_sp)
        real(sp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        real(sp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_sp
        do i = 1,m
           eye(i,i) = 1.0_sp
        end do
        a = eye + sample_matrix_s(m,m,shift=0.0_sp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_s

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: rect(2,3),z(3,3),empty(0,0)
        real(sp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_sp
        z = 0.0_sp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_s

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_s(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(sp),intent(in) :: shift
        real(sp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = sin(real(3*i + 5*j*j,kind=sp))
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_s

    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        real(dp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_dp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = i + j
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_d

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_dp
        d(1,1) = 1.0_dp
        d(2,2) = 1.0e-2_dp
        d(3,3) = 1.0e-4_dp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_dp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_dp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_dp,atol=1.0e-1_dp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_dp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_dp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_d

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(dp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_d(n,n,shift=1.0_dp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_dp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_dp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_dp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(dp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(dp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              real(dp) :: d(3,3)
              d = 0.0_dp
              d(1,1) = 1.0_dp
              d(2,2) = -2.0_dp
              d(3,3) = 4.0_dp
           end function diagonal

    end subroutine test_cond_d

    !> Inputs that cond rejects
    subroutine test_cond_invalid_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: a(3,3),rect(3,2),empty(0,0)
        real(dp) :: c
        type(la_state) :: state

        a = 1.0_dp
        rect = 1.0_dp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_d

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp),parameter :: tol_hilbert = 1.0e5_dp*epsilon(0.0_dp)
        real(dp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(dp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_d(n,n,shift=2.0_dp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_dp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(dp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=dp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_d

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(dp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative entry
        d3 = 0.0_dp
        d3(1,1) = 1.0_dp
        d3(2,2) = -2.0_dp
        d3(3,3) = 3.0_dp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad - log(6.0_dp)) <= tol,'negative determinant')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_dp
        do i = 1,50
           big(i,i) = 10.0_dp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_dp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_dp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_dp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_d

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(dp),parameter :: tol = 1000*epsilon(0.0_dp)
        real(dp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        real(dp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_dp
        do i = 1,m
           eye(i,i) = 1.0_dp
        end do
        a = eye + sample_matrix_d(m,m,shift=0.0_dp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_d

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: rect(2,3),z(3,3),empty(0,0)
        real(dp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_dp
        z = 0.0_dp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_d

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_d(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(dp),intent(in) :: shift
        real(dp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = sin(real(3*i + 5*j*j,kind=dp))
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_d

#ifdef LA_WITH_XDP
    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        real(xdp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_xdp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = i + j
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_x

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_xdp
        d(1,1) = 1.0_xdp
        d(2,2) = 1.0e-2_xdp
        d(3,3) = 1.0e-4_xdp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_xdp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_xdp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_xdp,atol=1.0e-1_xdp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_xdp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_xdp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_x

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        real(xdp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(xdp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_x(n,n,shift=1.0_xdp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_xdp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_xdp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_xdp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(xdp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(xdp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              real(xdp) :: d(3,3)
              d = 0.0_xdp
              d(1,1) = 1.0_xdp
              d(2,2) = -2.0_xdp
              d(3,3) = 4.0_xdp
           end function diagonal

    end subroutine test_cond_x

    !> Inputs that cond rejects
    subroutine test_cond_invalid_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp) :: a(3,3),rect(3,2),empty(0,0)
        real(xdp) :: c
        type(la_state) :: state

        a = 1.0_xdp
        rect = 1.0_xdp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_x

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        real(xdp),parameter :: tol_hilbert = 1.0e5_xdp*epsilon(0.0_xdp)
        real(xdp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(xdp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_x(n,n,shift=2.0_xdp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_xdp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(xdp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=xdp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_x

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        real(xdp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(xdp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative entry
        d3 = 0.0_xdp
        d3(1,1) = 1.0_xdp
        d3(2,2) = -2.0_xdp
        d3(3,3) = 3.0_xdp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad - log(6.0_xdp)) <= tol,'negative determinant')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_xdp
        do i = 1,50
           big(i,i) = 10.0_xdp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_xdp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_xdp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_xdp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_x

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(xdp),parameter :: tol = 1000*epsilon(0.0_xdp)
        real(xdp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        real(xdp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_xdp
        do i = 1,m
           eye(i,i) = 1.0_xdp
        end do
        a = eye + sample_matrix_x(m,m,shift=0.0_xdp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_x

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp) :: rect(2,3),z(3,3),empty(0,0)
        real(xdp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_xdp
        z = 0.0_xdp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_x

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_x(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(xdp),intent(in) :: shift
        real(xdp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = sin(real(3*i + 5*j*j,kind=xdp))
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_x
#endif

#ifdef LA_WITH_QP
    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        real(qp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_qp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = i + j
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_q

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_qp
        d(1,1) = 1.0_qp
        d(2,2) = 1.0e-2_qp
        d(3,3) = 1.0e-4_qp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_qp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_qp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_qp,atol=1.0e-1_qp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_qp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_qp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_q

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(qp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_q(n,n,shift=1.0_qp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_qp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_qp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_qp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(qp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(qp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              real(qp) :: d(3,3)
              d = 0.0_qp
              d(1,1) = 1.0_qp
              d(2,2) = -2.0_qp
              d(3,3) = 4.0_qp
           end function diagonal

    end subroutine test_cond_q

    !> Inputs that cond rejects
    subroutine test_cond_invalid_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: a(3,3),rect(3,2),empty(0,0)
        real(qp) :: c
        type(la_state) :: state

        a = 1.0_qp
        rect = 1.0_qp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_q

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp),parameter :: tol_hilbert = 1.0e5_qp*epsilon(0.0_qp)
        real(qp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(qp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_q(n,n,shift=2.0_qp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_qp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(qp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=qp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_q

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(qp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative entry
        d3 = 0.0_qp
        d3(1,1) = 1.0_qp
        d3(2,2) = -2.0_qp
        d3(3,3) = 3.0_qp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad - log(6.0_qp)) <= tol,'negative determinant')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_qp
        do i = 1,50
           big(i,i) = 10.0_qp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_qp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_qp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_qp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_q

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(qp),parameter :: tol = 1000*epsilon(0.0_qp)
        real(qp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        real(qp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_qp
        do i = 1,m
           eye(i,i) = 1.0_qp
        end do
        a = eye + sample_matrix_q(m,m,shift=0.0_qp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_q

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: rect(2,3),z(3,3),empty(0,0)
        real(qp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_qp
        z = 0.0_qp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_q

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_q(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(qp),intent(in) :: shift
        real(qp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = sin(real(3*i + 5*j*j,kind=qp))
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_q
#endif

    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        complex(sp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_sp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = cmplx(i + j,j - i,kind=sp)
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_c

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_c(error)
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_sp
        d(1,1) = 1.0_sp
        d(2,2) = 1.0e-2_sp
        d(3,3) = 1.0e-4_sp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_sp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_sp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_sp,atol=1.0e-1_sp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_sp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_sp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_c

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(sp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_c(n,n,shift=1.0_sp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_sp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_sp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_sp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(sp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(sp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              complex(sp) :: d(3,3)
              d = 0.0_sp
              d(1,1) = 1.0_sp
              d(2,2) = -2.0_sp
              d(3,3) = 4.0_sp
           end function diagonal

    end subroutine test_cond_c

    !> Inputs that cond rejects
    subroutine test_cond_invalid_c(error)
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: a(3,3),rect(3,2),empty(0,0)
        real(sp) :: c
        type(la_state) :: state

        a = 1.0_sp
        rect = 1.0_sp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_c

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp),parameter :: tol_hilbert = 1.0e5_sp*epsilon(0.0_sp)
        complex(sp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(sp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_c(n,n,shift=2.0_sp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_sp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(sp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=sp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_c

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        complex(sp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(sp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative and an imaginary entry
        d3 = 0.0_sp
        d3(1,1) = (0.0_sp,1.0_sp)
        d3(2,2) = 2.0_sp
        d3(3,3) = -1.0_sp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn - (0.0_sp,-1.0_sp)) <= tol .and. abs(lad - log(2.0_sp)) <= tol, &
                   'complex sign')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_sp
        do i = 1,50
           big(i,i) = 10.0_sp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_sp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_sp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_sp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_c

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(sp),parameter :: tol = 1000*epsilon(0.0_sp)
        complex(sp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        complex(sp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_sp
        do i = 1,m
           eye(i,i) = 1.0_sp
        end do
        a = eye + sample_matrix_c(m,m,shift=0.0_sp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_c

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_c(error)
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: rect(2,3),z(3,3),empty(0,0)
        complex(sp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_sp
        z = 0.0_sp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_c

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_c(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(sp),intent(in) :: shift
        complex(sp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = cmplx(sin(real(3*i + 5*j*j,kind=sp)),cos(real(2*i*i + 7*j,kind=sp)),kind=sp)
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_c

    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        complex(dp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_dp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = cmplx(i + j,j - i,kind=dp)
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_z

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_z(error)
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_dp
        d(1,1) = 1.0_dp
        d(2,2) = 1.0e-2_dp
        d(3,3) = 1.0e-4_dp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_dp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_dp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_dp,atol=1.0e-1_dp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_dp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_dp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_z

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(dp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_z(n,n,shift=1.0_dp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_dp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_dp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_dp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(dp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(dp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              complex(dp) :: d(3,3)
              d = 0.0_dp
              d(1,1) = 1.0_dp
              d(2,2) = -2.0_dp
              d(3,3) = 4.0_dp
           end function diagonal

    end subroutine test_cond_z

    !> Inputs that cond rejects
    subroutine test_cond_invalid_z(error)
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: a(3,3),rect(3,2),empty(0,0)
        real(dp) :: c
        type(la_state) :: state

        a = 1.0_dp
        rect = 1.0_dp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_z

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp),parameter :: tol_hilbert = 1.0e5_dp*epsilon(0.0_dp)
        complex(dp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(dp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_z(n,n,shift=2.0_dp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_dp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(dp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=dp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_z

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        complex(dp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(dp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative and an imaginary entry
        d3 = 0.0_dp
        d3(1,1) = (0.0_dp,1.0_dp)
        d3(2,2) = 2.0_dp
        d3(3,3) = -1.0_dp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn - (0.0_dp,-1.0_dp)) <= tol .and. abs(lad - log(2.0_dp)) <= tol, &
                   'complex sign')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_dp
        do i = 1,50
           big(i,i) = 10.0_dp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_dp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_dp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_dp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_z

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(dp),parameter :: tol = 1000*epsilon(0.0_dp)
        complex(dp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        complex(dp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_dp
        do i = 1,m
           eye(i,i) = 1.0_dp
        end do
        a = eye + sample_matrix_z(m,m,shift=0.0_dp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_z

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_z(error)
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: rect(2,3),z(3,3),empty(0,0)
        complex(dp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_dp
        z = 0.0_dp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_z

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_z(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(dp),intent(in) :: shift
        complex(dp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = cmplx(sin(real(3*i + 5*j*j,kind=dp)),cos(real(2*i*i + 7*j,kind=dp)),kind=dp)
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_z

#ifdef LA_WITH_XDP
    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        complex(xdp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_xdp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = cmplx(i + j,j - i,kind=xdp)
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_y

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_y(error)
        type(error_type),allocatable,intent(out) :: error

        complex(xdp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_xdp
        d(1,1) = 1.0_xdp
        d(2,2) = 1.0e-2_xdp
        d(3,3) = 1.0e-4_xdp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_xdp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_xdp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_xdp,atol=1.0e-1_xdp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_xdp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_xdp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_y

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        complex(xdp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(xdp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_y(n,n,shift=1.0_xdp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_xdp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_xdp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_xdp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(xdp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(xdp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              complex(xdp) :: d(3,3)
              d = 0.0_xdp
              d(1,1) = 1.0_xdp
              d(2,2) = -2.0_xdp
              d(3,3) = 4.0_xdp
           end function diagonal

    end subroutine test_cond_y

    !> Inputs that cond rejects
    subroutine test_cond_invalid_y(error)
        type(error_type),allocatable,intent(out) :: error

        complex(xdp) :: a(3,3),rect(3,2),empty(0,0)
        real(xdp) :: c
        type(la_state) :: state

        a = 1.0_xdp
        rect = 1.0_xdp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_y

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        real(xdp),parameter :: tol_hilbert = 1.0e5_xdp*epsilon(0.0_xdp)
        complex(xdp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(xdp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_y(n,n,shift=2.0_xdp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_xdp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(xdp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=xdp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_y

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
        complex(xdp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(xdp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative and an imaginary entry
        d3 = 0.0_xdp
        d3(1,1) = (0.0_xdp,1.0_xdp)
        d3(2,2) = 2.0_xdp
        d3(3,3) = -1.0_xdp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn - (0.0_xdp,-1.0_xdp)) <= tol .and. abs(lad - log(2.0_xdp)) <= tol, &
                   'complex sign')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_xdp
        do i = 1,50
           big(i,i) = 10.0_xdp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_xdp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_xdp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_xdp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_y

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(xdp),parameter :: tol = 1000*epsilon(0.0_xdp)
        complex(xdp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        complex(xdp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_xdp
        do i = 1,m
           eye(i,i) = 1.0_xdp
        end do
        a = eye + sample_matrix_y(m,m,shift=0.0_xdp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_y

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_y(error)
        type(error_type),allocatable,intent(out) :: error

        complex(xdp) :: rect(2,3),z(3,3),empty(0,0)
        complex(xdp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_xdp
        z = 0.0_xdp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_y

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_y(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(xdp),intent(in) :: shift
        complex(xdp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = cmplx(sin(real(3*i + 5*j*j,kind=xdp)),cos(real(2*i*i + 7*j,kind=xdp)),kind=xdp)
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_y
#endif

#ifdef LA_WITH_QP
    !> Products of rank r = 0, ..., n of integer factors, so the products are exact: the first r
    !> columns of a Vandermonde matrix with distinct nodes, by the first r rows of an upper
    !> triangular matrix with a nonzero diagonal
    subroutine test_matrix_rank_ladder_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 6_ilp
        integer(ilp),parameter :: n = 4_ilp
        complex(qp) :: u(m,n),v(n,n),a(m,n)
        integer(ilp) :: r,i,j

        do concurrent(i=1:m,j=1:n)
           u(i,j) = i**(j - 1)
        end do
        v = 0.0_qp
        do concurrent(i=1:n,j=1:n,j >= i)
           v(i,j) = cmplx(i + j,j - i,kind=qp)
        end do

        do r = 0,n
           a = matmul(u(:,1:r),v(1:r,:))

           call check(error,matrix_rank(a) == r,'rank of a tall matrix')
           if (allocated(error)) return

           call check(error,matrix_rank(transpose(a)) == r,'rank of a wide matrix')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_rank_ladder_w

    !> diag(1, 1e-2, 1e-4) with different tolerances, an empty matrix, invalid tolerances
    subroutine test_matrix_rank_tolerances_w(error)
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: d(3,3),empty(0,3)
        integer(ilp) :: r
        type(la_state) :: state

        d = 0.0_qp
        d(1,1) = 1.0_qp
        d(2,2) = 1.0e-2_qp
        d(3,3) = 1.0e-4_qp

        call check(error,matrix_rank(d) == 3,'default tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=1.0e-3_qp) == 2,'absolute tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-1_qp) == 1,'relative tolerance')
        if (allocated(error)) return

        call check(error,matrix_rank(d,rtol=1.0e-3_qp,atol=1.0e-1_qp) == 1,'larger of the two tolerances')
        if (allocated(error)) return

        call check(error,matrix_rank(d,atol=2.0_qp) == 0,'tolerance above every singular value')
        if (allocated(error)) return

        r = matrix_rank(empty,err=state)
        call check(error,state%ok() .and. r == 0,'empty matrix')
        if (allocated(error)) return

        r = matrix_rank(d,rtol=-1.0_qp,err=state)
        call check(error,state%error() .and. r == 0,'negative tolerance')
        if (allocated(error)) return

    end subroutine test_matrix_rank_tolerances_w

    !> Condition numbers against the norms of the matrix and of its inverse
    subroutine test_cond_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5_ilp
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),ainv(n,n),rect(n,n - 2)
        real(qp) :: s(n),c,ref
        type(la_state) :: state

        a = sample_matrix_w(n,n,shift=1.0_qp)
        ainv = inv(a)

        ! 2-norm: ratio of the extreme singular values
        s = svdvals(a)
        ref = s(1)/s(n)

        call check(error,is_close(cond(a),ref),'default (2-norm) condition number')
        if (allocated(error)) return

        call check(error,is_close(cond(a,2_ilp),ref) .and. is_close(cond(a,'2'),ref),'order=2')
        if (allocated(error)) return

        c = cond(a,err=state)
        call check(error,state%ok() .and. is_close(c,ref),'default condition number with err')
        if (allocated(error)) return

        ! 1-norm and infinity norm: gecon estimates the norm of the inverse from below, usually
        ! within a factor of 3
        ref = mnorm(a,'1')*mnorm(ainv,'1')
        call check(error,is_estimate(cond(a,1_ilp),ref) .and. cond(a,'1') == cond(a,1_ilp),'1-norm')
        if (allocated(error)) return

        ref = mnorm(a,'inf')*mnorm(ainv,'inf')
        call check(error,is_estimate(cond(a,huge(0_ilp)),ref) .and. cond(a,'inf') == cond(a,huge(0_ilp)), &
                   'infinity norm')
        if (allocated(error)) return

        c = cond(a,'inf',err=state)
        call check(error,state%ok() .and. is_estimate(c,ref),'infinity norm with err')
        if (allocated(error)) return

        ! A diagonal matrix, whose estimate is exact
        call check(error,is_close(cond(diagonal(),1_ilp),4.0_qp) .and. &
                          is_close(cond(diagonal(),'inf'),4.0_qp),'diagonal matrix')
        if (allocated(error)) return

        ! Frobenius norm
        ref = mnorm(a,'fro')*mnorm(ainv,'fro')
        call check(error,is_close(cond(a,'fro'),ref),'Frobenius norm')
        if (allocated(error)) return

        ! Rectangular matrix, 2-norm
        rect = a(:,1:n - 2)
        s(1:n - 2) = svdvals(rect)
        call check(error,is_close(cond(rect),s(1)/s(n - 2)),'rectangular matrix')
        if (allocated(error)) return

        ! A zero column: infinite condition number, not an error
        a(:,2) = 0.0_qp
        c = cond(a,1_ilp,err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(a,'fro',err=state)
        call check(error,state%ok() .and. .not. ieee_is_finite(c) .and. c > 0,'singular matrix, Frobenius norm')
        if (allocated(error)) return

        contains

           elemental logical function is_close(x,y)
              real(qp),intent(in) :: x,y
              is_close = abs(x - y) <= tol*abs(y)
           end function is_close

           elemental logical function is_estimate(x,y)
              real(qp),intent(in) :: x,y
              is_estimate = x <= (1 + tol)*y .and. x >= y/3
           end function is_estimate

           pure function diagonal() result(d)
              complex(qp) :: d(3,3)
              d = 0.0_qp
              d(1,1) = 1.0_qp
              d(2,2) = -2.0_qp
              d(3,3) = 4.0_qp
           end function diagonal

    end subroutine test_cond_w

    !> Inputs that cond rejects
    subroutine test_cond_invalid_w(error)
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: a(3,3),rect(3,2),empty(0,0)
        real(qp) :: c
        type(la_state) :: state

        a = 1.0_qp
        rect = 1.0_qp

        c = cond(rect,1_ilp,err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, 1-norm')
        if (allocated(error)) return

        c = cond(rect,'fro',err=state)
        call check(error,state%error() .and. c == 0,'rectangular matrix, Frobenius norm')
        if (allocated(error)) return

        c = cond(a,3_ilp,err=state)
        call check(error,state%error() .and. c == 0,'order=3')
        if (allocated(error)) return

        c = cond(a,'-inf',err=state)
        call check(error,state%error() .and. c == 0,'order=-inf')
        if (allocated(error)) return

        c = cond(a,'nuc',err=state)
        call check(error,state%error() .and. c == 0,'order=nuc')
        if (allocated(error)) return

        c = cond(empty,err=state)
        call check(error,state%error() .and. c == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_cond_invalid_w

    !> Log-determinant against det, and of the Hilbert matrix against its closed form
    subroutine test_slogdet_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 4_ilp
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp),parameter :: tol_hilbert = 1.0e5_qp*epsilon(0.0_qp)
        complex(qp) :: a(n,n),aorig(n,n),h(n,n),sgn,sgn2,d
        real(qp) :: lad,lad2,ref
        integer(ilp) :: i,j
        type(la_state) :: state

        ! A well-scaled matrix
        a = sample_matrix_w(n,n,shift=2.0_qp*n)
        aorig = a
        d = det(a)

        call slogdet(a,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn*exp(lad) - d) <= tol*abs(d),'sign*exp(logabsdet) = det(a)')
        if (allocated(error)) return

        call check(error,abs(abs(sgn) - 1) <= tol,'unit modulus sign')
        if (allocated(error)) return

        call check(error,all(a == aorig),'a is not overwritten by default')
        if (allocated(error)) return

        call slogdet(a,sgn2,lad2,overwrite_a=.true.,err=state)
        call check(error,state%ok() .and. abs(sgn2 - sgn) <= tol .and. abs(lad2 - lad) <= tol*abs(lad), &
                   'same result with overwrite_a')
        if (allocated(error)) return

        ! Hilbert matrix: det(H) = c(n)**4/c(2n), with c(k) = 1! 2! ... (k-1)!
        do concurrent(i=1:n,j=1:n)
           h(i,j) = 1.0_qp/(i + j - 1)
        end do
        ref = 4*log_superfactorial(n) - log_superfactorial(2*n)

        call slogdet(h,sgn,lad,err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        call check(error,abs(sgn - 1) <= tol .and. abs(lad - ref) <= tol_hilbert*abs(ref),'Hilbert matrix')
        if (allocated(error)) return

        contains

           !> log(1! 2! ... (k-1)!)
           pure real(qp) function log_superfactorial(k)
              integer(ilp),intent(in) :: k
              integer(ilp) :: i
              log_superfactorial = sum([(log_gamma(real(i + 1,kind=qp)),i=1,k - 1)])
           end function log_superfactorial

    end subroutine test_slogdet_w

    !> Signs, a determinant out of range, singular, empty and non-square matrices
    subroutine test_slogdet_special_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        complex(qp) :: p(2,2),d3(3,3),big(50,50),z(3,3),empty(0,0),rect(2,3),sgn
        real(qp) :: lad
        integer(ilp) :: i
        type(la_state) :: state

        ! A row interchange: det = -1
        p = reshape([0,1,1,0], [2,2])
        call slogdet(p,sgn,lad)
        call check(error,abs(sgn + 1) <= tol .and. abs(lad) <= tol,'permutation matrix')
        if (allocated(error)) return

        ! Diagonal matrix with a negative and an imaginary entry
        d3 = 0.0_qp
        d3(1,1) = (0.0_qp,1.0_qp)
        d3(2,2) = 2.0_qp
        d3(3,3) = -1.0_qp
        call slogdet(d3,sgn,lad)
        call check(error,abs(sgn - (0.0_qp,-1.0_qp)) <= tol .and. abs(lad - log(2.0_qp)) <= tol, &
                   'complex sign')
        if (allocated(error)) return

        ! det = 10**50, beyond the range of real(sp)
        big = 0.0_qp
        do i = 1,50
           big(i,i) = 10.0_qp
        end do
        call slogdet(big,sgn,lad)
        call check(error,abs(sgn - 1) <= tol .and. abs(lad - 50*log(10.0_qp)) <= tol*lad,'determinant out of range')
        if (allocated(error)) return

        ! Singular matrix: sign 0, logabsdet -inf, no error
        z = 0.0_qp
        call slogdet(z,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 0 .and. .not. ieee_is_finite(lad) .and. ieee_is_negative(lad), &
                   'singular matrix')
        if (allocated(error)) return

        ! Empty matrix: det = 1
        call slogdet(empty,sgn,lad,err=state)
        call check(error,state%ok() .and. sgn == 1 .and. lad == 0,'empty matrix')
        if (allocated(error)) return

        ! Non-square matrix
        rect = 1.0_qp
        call slogdet(rect,sgn,lad,err=state)
        call check(error,state%error(),'non-square matrix')
        if (allocated(error)) return

    end subroutine test_slogdet_special_w

    !> Matrix powers against repeated products, for positive, zero and negative powers
    subroutine test_matrix_power_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: m = 4_ilp
        real(qp),parameter :: tol = 1000*epsilon(0.0_qp)
        complex(qp) :: a(m,m),ainv(m,m),eye(m,m),ref(m,m)
        complex(qp),allocatable :: an(:,:)
        integer(ilp) :: i,k
        type(la_state) :: state

        ! A well-conditioned matrix: the identity plus a perturbation of norm below 1
        eye = 0.0_qp
        do i = 1,m
           eye(i,i) = 1.0_qp
        end do
        a = eye + sample_matrix_w(m,m,shift=0.0_qp)/(2*m)
        ainv = inv(a)

        ! Zero power: the identity
        an = matrix_power(a,0_ilp)
        call check(error,all(an == eye),'a**0')
        if (allocated(error)) return

        ! Positive powers
        ref = eye
        do k = 1,7
           ref = matmul(ref,a)
           an = matrix_power(a,k)
           call check(error,all(abs(an - ref) <= tol*maxval(abs(ref))),'positive power')
           if (allocated(error)) return
        end do

        ! Negative powers
        ref = eye
        do k = 1,3
           ref = matmul(ref,ainv)
           an = matrix_power(a,-k,err=state)
           call check(error,state%ok() .and. all(abs(an - ref) <= tol*maxval(abs(ref))),'negative power')
           if (allocated(error)) return
        end do

    end subroutine test_matrix_power_w

    !> Inputs that matrix_power rejects or treats specially
    subroutine test_matrix_power_invalid_w(error)
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: rect(2,3),z(3,3),empty(0,0)
        complex(qp),allocatable :: an(:,:)
        type(la_state) :: state

        rect = 1.0_qp
        z = 0.0_qp

        an = matrix_power(rect,2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'non-square matrix')
        if (allocated(error)) return

        an = matrix_power(z,-2_ilp,err=state)
        call check(error,state%error() .and. size(an) == 0,'negative power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,0_ilp,err=state)
        call check(error,state%ok() .and. all(shape(an) == [3,3]),'zero power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(z,3_ilp,err=state)
        call check(error,state%ok() .and. all(an == 0),'positive power of a singular matrix')
        if (allocated(error)) return

        an = matrix_power(empty,-3_ilp,err=state)
        call check(error,state%ok() .and. size(an) == 0,'empty matrix')
        if (allocated(error)) return

    end subroutine test_matrix_power_invalid_w

    !> A deterministic matrix with entries in [-1,1] (each part, if complex), plus shift*I
    pure function sample_matrix_w(m,n,shift) result(a)
        integer(ilp),intent(in) :: m,n
        real(qp),intent(in) :: shift
        complex(qp) :: a(m,n)

        integer(ilp) :: i,j

        do concurrent(i=1:m,j=1:n)
           a(i,j) = cmplx(sin(real(3*i + 5*j*j,kind=qp)),cos(real(2*i*i + 7*j,kind=qp)),kind=qp)
           if (i == j) a(i,j) = a(i,j) + shift
        end do

    end function sample_matrix_w
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

end module test_linalg_scalars

program test_scalars
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_scalars,only:test_matrix_scalars
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_scalars",test_matrix_scalars) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_scalars
