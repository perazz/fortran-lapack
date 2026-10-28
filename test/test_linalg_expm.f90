! Test the matrix exponential
module test_linalg_expm
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:expm,matrix_exp,eye,norm
    use la_state_type,only:la_state

    implicit none(type,external)
    private

    public :: test_expm_computation

    contains

    !> Matrix exponential tests
    subroutine test_expm_computation(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("expm_s",test_expm_s))
        call add_test(tests,new_unittest("error_handling_expm_s",test_error_handling_expm_s))
        call add_test(tests,new_unittest("expm_rotation_s",test_expm_rotation_s))
        call add_test(tests,new_unittest("expm_err_form_s",test_expm_err_form_s))
        call add_test(tests,new_unittest("expm_pure_s",test_expm_pure_s))
        call add_test(tests,new_unittest("expm_zero_and_empty_s",test_expm_zero_and_empty_s))
        call add_test(tests,new_unittest("expm_d",test_expm_d))
        call add_test(tests,new_unittest("error_handling_expm_d",test_error_handling_expm_d))
        call add_test(tests,new_unittest("expm_rotation_d",test_expm_rotation_d))
        call add_test(tests,new_unittest("expm_err_form_d",test_expm_err_form_d))
        call add_test(tests,new_unittest("expm_pure_d",test_expm_pure_d))
        call add_test(tests,new_unittest("expm_zero_and_empty_d",test_expm_zero_and_empty_d))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("expm_x",test_expm_x))
        call add_test(tests,new_unittest("error_handling_expm_x",test_error_handling_expm_x))
        call add_test(tests,new_unittest("expm_rotation_x",test_expm_rotation_x))
        call add_test(tests,new_unittest("expm_err_form_x",test_expm_err_form_x))
        call add_test(tests,new_unittest("expm_pure_x",test_expm_pure_x))
        call add_test(tests,new_unittest("expm_zero_and_empty_x",test_expm_zero_and_empty_x))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("expm_q",test_expm_q))
        call add_test(tests,new_unittest("error_handling_expm_q",test_error_handling_expm_q))
        call add_test(tests,new_unittest("expm_rotation_q",test_expm_rotation_q))
        call add_test(tests,new_unittest("expm_err_form_q",test_expm_err_form_q))
        call add_test(tests,new_unittest("expm_pure_q",test_expm_pure_q))
        call add_test(tests,new_unittest("expm_zero_and_empty_q",test_expm_zero_and_empty_q))
#endif
        call add_test(tests,new_unittest("expm_c",test_expm_c))
        call add_test(tests,new_unittest("error_handling_expm_c",test_error_handling_expm_c))
        call add_test(tests,new_unittest("expm_rotation_c",test_expm_rotation_c))
        call add_test(tests,new_unittest("expm_err_form_c",test_expm_err_form_c))
        call add_test(tests,new_unittest("expm_pure_c",test_expm_pure_c))
        call add_test(tests,new_unittest("expm_zero_and_empty_c",test_expm_zero_and_empty_c))
        call add_test(tests,new_unittest("expm_z",test_expm_z))
        call add_test(tests,new_unittest("error_handling_expm_z",test_error_handling_expm_z))
        call add_test(tests,new_unittest("expm_rotation_z",test_expm_rotation_z))
        call add_test(tests,new_unittest("expm_err_form_z",test_expm_err_form_z))
        call add_test(tests,new_unittest("expm_pure_z",test_expm_pure_z))
        call add_test(tests,new_unittest("expm_zero_and_empty_z",test_expm_zero_and_empty_z))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("expm_y",test_expm_y))
        call add_test(tests,new_unittest("error_handling_expm_y",test_error_handling_expm_y))
        call add_test(tests,new_unittest("expm_rotation_y",test_expm_rotation_y))
        call add_test(tests,new_unittest("expm_err_form_y",test_expm_err_form_y))
        call add_test(tests,new_unittest("expm_pure_y",test_expm_pure_y))
        call add_test(tests,new_unittest("expm_zero_and_empty_y",test_expm_zero_and_empty_y))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("expm_w",test_expm_w))
        call add_test(tests,new_unittest("error_handling_expm_w",test_error_handling_expm_w))
        call add_test(tests,new_unittest("expm_rotation_w",test_expm_rotation_w))
        call add_test(tests,new_unittest("expm_err_form_w",test_expm_err_form_w))
        call add_test(tests,new_unittest("expm_pure_w",test_expm_pure_w))
        call add_test(tests,new_unittest("expm_zero_and_empty_w",test_expm_zero_and_empty_w))
#endif

    end subroutine test_expm_computation

    !> Matrix exponential with analytic expression
    subroutine test_expm_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        real(sp) :: A(n,n),E(n,n),Eref(n,n)
        real(sp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_sp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_sp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_sp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_sp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_s

    !> Test error handler
    subroutine test_error_handling_expm_s(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        real(sp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_sp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_sp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_s

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: t = 3.0_sp
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_sp,t,-t,0.0_sp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_s

    !> The function form with a state flag
    subroutine test_expm_err_form_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: t = 3.0_sp
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp) :: A(2,2),Eref(2,2),B(2,3)
        real(sp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_sp,t,-t,0.0_sp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_sp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_s

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_sp,1.0_sp,-1.0_sp,0.0_sp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_s(real(k,sp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_s(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_s

    pure function squared_exp_s(A) result(E2)
        real(sp),intent(in) :: A(:,:)
        real(sp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_s

    pure function rotation_s(t) result(R)
        real(sp),intent(in) :: t
        real(sp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_s

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_s(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        real(sp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        real(sp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_sp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_sp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_s

    !> Matrix exponential with analytic expression
    subroutine test_expm_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        real(dp) :: A(n,n),E(n,n),Eref(n,n)
        real(dp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_dp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_dp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_dp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_dp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_d

    !> Test error handler
    subroutine test_error_handling_expm_d(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        real(dp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_dp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_dp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_d

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: t = 3.0_dp
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_dp,t,-t,0.0_dp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_d

    !> The function form with a state flag
    subroutine test_expm_err_form_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: t = 3.0_dp
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp) :: A(2,2),Eref(2,2),B(2,3)
        real(dp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_dp,t,-t,0.0_dp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_dp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_d

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_d(real(k,dp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_d(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_d

    pure function squared_exp_d(A) result(E2)
        real(dp),intent(in) :: A(:,:)
        real(dp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_d

    pure function rotation_d(t) result(R)
        real(dp),intent(in) :: t
        real(dp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_d

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_d(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        real(dp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        real(dp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_dp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_dp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_d

#ifdef LA_WITH_XDP
    !> Matrix exponential with analytic expression
    subroutine test_expm_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        real(xdp) :: A(n,n),E(n,n),Eref(n,n)
        real(xdp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_xdp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_xdp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_xdp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_xdp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_x

    !> Test error handler
    subroutine test_error_handling_expm_x(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        real(xdp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_xdp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_xdp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_x

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: t = 3.0_xdp
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_xdp,t,-t,0.0_xdp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_x

    !> The function form with a state flag
    subroutine test_expm_err_form_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: t = 3.0_xdp
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp) :: A(2,2),Eref(2,2),B(2,3)
        real(xdp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_xdp,t,-t,0.0_xdp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_xdp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_x

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_xdp,1.0_xdp,-1.0_xdp,0.0_xdp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_x(real(k,xdp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_x(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_x

    pure function squared_exp_x(A) result(E2)
        real(xdp),intent(in) :: A(:,:)
        real(xdp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_x

    pure function rotation_x(t) result(R)
        real(xdp),intent(in) :: t
        real(xdp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_x

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_x(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        real(xdp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        real(xdp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_xdp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_xdp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_x
#endif

#ifdef LA_WITH_QP
    !> Matrix exponential with analytic expression
    subroutine test_expm_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        real(qp) :: A(n,n),E(n,n),Eref(n,n)
        real(qp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_qp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_qp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_qp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_qp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_q

    !> Test error handler
    subroutine test_error_handling_expm_q(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        real(qp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_qp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_qp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_q

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: t = 3.0_qp
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_qp,t,-t,0.0_qp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_q

    !> The function form with a state flag
    subroutine test_expm_err_form_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: t = 3.0_qp
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp) :: A(2,2),Eref(2,2),B(2,3)
        real(qp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_qp,t,-t,0.0_qp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_qp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_q

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_qp,1.0_qp,-1.0_qp,0.0_qp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_q(real(k,qp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_q(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_q

    pure function squared_exp_q(A) result(E2)
        real(qp),intent(in) :: A(:,:)
        real(qp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_q

    pure function rotation_q(t) result(R)
        real(qp),intent(in) :: t
        real(qp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_q

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_q(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        real(qp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        real(qp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_qp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_qp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_q
#endif

    !> Matrix exponential with analytic expression
    subroutine test_expm_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        complex(sp) :: A(n,n),E(n,n),Eref(n,n)
        real(sp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_sp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_sp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_sp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_sp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_c

    !> Test error handler
    subroutine test_error_handling_expm_c(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        complex(sp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_sp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_sp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_c

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: t = 3.0_sp
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        complex(sp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_sp,t,-t,0.0_sp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_c

    !> The function form with a state flag
    subroutine test_expm_err_form_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: t = 3.0_sp
        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        complex(sp) :: A(2,2),Eref(2,2),B(2,3)
        complex(sp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_sp,t,-t,0.0_sp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_sp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_c

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        complex(sp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_sp,1.0_sp,-1.0_sp,0.0_sp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_c(real(k,sp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_c(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_c

    pure function squared_exp_c(A) result(E2)
        complex(sp),intent(in) :: A(:,:)
        complex(sp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_c

    pure function rotation_c(t) result(R)
        real(sp),intent(in) :: t
        complex(sp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_c

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_c(error)
        type(error_type),allocatable,intent(out) :: error

        real(sp),parameter :: tol = 1000*epsilon(1.0_sp)
        complex(sp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        complex(sp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_sp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_sp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_c

    !> Matrix exponential with analytic expression
    subroutine test_expm_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        complex(dp) :: A(n,n),E(n,n),Eref(n,n)
        real(dp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_dp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_dp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_dp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_dp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_z

    !> Test error handler
    subroutine test_error_handling_expm_z(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        complex(dp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_dp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_dp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_z

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: t = 3.0_dp
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        complex(dp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_dp,t,-t,0.0_dp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_z

    !> The function form with a state flag
    subroutine test_expm_err_form_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: t = 3.0_dp
        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        complex(dp) :: A(2,2),Eref(2,2),B(2,3)
        complex(dp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_dp,t,-t,0.0_dp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_dp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_z

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        complex(dp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_z(real(k,dp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_z(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_z

    pure function squared_exp_z(A) result(E2)
        complex(dp),intent(in) :: A(:,:)
        complex(dp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_z

    pure function rotation_z(t) result(R)
        real(dp),intent(in) :: t
        complex(dp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_z

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_z(error)
        type(error_type),allocatable,intent(out) :: error

        real(dp),parameter :: tol = 1000*epsilon(1.0_dp)
        complex(dp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        complex(dp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_dp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_dp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_z

#ifdef LA_WITH_XDP
    !> Matrix exponential with analytic expression
    subroutine test_expm_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        complex(xdp) :: A(n,n),E(n,n),Eref(n,n)
        real(xdp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_xdp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_xdp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_xdp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_xdp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_y

    !> Test error handler
    subroutine test_error_handling_expm_y(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        complex(xdp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_xdp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_xdp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_y

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: t = 3.0_xdp
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        complex(xdp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_xdp,t,-t,0.0_xdp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_y

    !> The function form with a state flag
    subroutine test_expm_err_form_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: t = 3.0_xdp
        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        complex(xdp) :: A(2,2),Eref(2,2),B(2,3)
        complex(xdp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_xdp,t,-t,0.0_xdp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_xdp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_y

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        complex(xdp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_xdp,1.0_xdp,-1.0_xdp,0.0_xdp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_y(real(k,xdp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_y(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_y

    pure function squared_exp_y(A) result(E2)
        complex(xdp),intent(in) :: A(:,:)
        complex(xdp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_y

    pure function rotation_y(t) result(R)
        real(xdp),intent(in) :: t
        complex(xdp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_y

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_y(error)
        type(error_type),allocatable,intent(out) :: error

        real(xdp),parameter :: tol = 1000*epsilon(1.0_xdp)
        complex(xdp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        complex(xdp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_xdp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_xdp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_y
#endif

#ifdef LA_WITH_QP
    !> Matrix exponential with analytic expression
    subroutine test_expm_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        integer(ilp) :: i,j
        complex(qp) :: A(n,n),E(n,n),Eref(n,n)
        real(qp) :: err

        ! Nilpotent matrix with m on the first superdiagonal
        A = 0.0_qp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_qp
        end do

        ! Reference with analytical exponential
        Eref = eye(n,mold=1.0_qp)
        do i = 1,n - 1
            do j = 1,n - i
                Eref(i,i + j) = Eref(i,i + j - 1)*m/j
            end do
        end do

        ! Compute matrix exponential
        E = expm(A)

        ! Check result
        err = norm(Eref - E,'inf')
        call check(error,err < (n**2)*epsilon(1.0_qp),'Analytical matrix exponential.')
        if (allocated(error)) return

    end subroutine test_expm_w

    !> Test error handler
    subroutine test_error_handling_expm_w(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 5,m = 6
        type(la_state) :: err
        integer(ilp) :: i
        complex(qp) :: A(n,n),E(n,n)

        ! Initialize matrix
        A = 0.0_qp
        do i = 1,n - 1
            A(i,i + 1) = m*1.0_qp
        end do

        call matrix_exp(A,E,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,order=-1,err=err)
        call check(error,err%error(),'Negative Pade order')
        if (allocated(error)) return

        call matrix_exp(A,E(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

        call matrix_exp(A(:n,:n - 1),err=err)
        call check(error,err%error(),'Invalid matrix size')
        if (allocated(error)) return

    end subroutine test_error_handling_expm_w

    !> exp of a rotation generator: large enough to need scaling and squaring
    subroutine test_expm_rotation_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: t = 3.0_qp
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        complex(qp) :: A(2,2),E(2,2),Eref(2,2),Ein(2,2)
        type(la_state) :: state

        A = reshape([0.0_qp,t,-t,0.0_qp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A)
        call check(error,maxval(abs(E - Eref)) < tol,'rotation: expm(A)')
        if (allocated(error)) return

        call matrix_exp(A,E,err=state)
        call check(error,state%ok() .and. maxval(abs(E - Eref)) < tol,'rotation: matrix_exp(A, E)')
        if (allocated(error)) return

        Ein = A
        call matrix_exp(Ein,err=state)
        call check(error,state%ok() .and. maxval(abs(Ein - Eref)) < tol,'rotation: matrix_exp(A) in place')
        if (allocated(error)) return

    end subroutine test_expm_rotation_w

    !> The function form with a state flag
    subroutine test_expm_err_form_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: t = 3.0_qp
        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        complex(qp) :: A(2,2),Eref(2,2),B(2,3)
        complex(qp),allocatable :: E(:,:)
        type(la_state) :: state

        A = reshape([0.0_qp,t,-t,0.0_qp], [2,2])
        Eref = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

        E = expm(A,err=state)
        call check(error,state%ok(),'expm(A, err): '//state%print())
        if (allocated(error)) return
        call check(error,maxval(abs(E - Eref)) < tol,'expm(A, err): result')
        if (allocated(error)) return

        ! The order reaches the computation: same result as the pure form of that order
        E = expm(A,order=6,err=state)
        call check(error,state%ok() .and. maxval(abs(E - expm(A,order=6))) < tol,'expm(A, order=6, err)')
        if (allocated(error)) return

        E = expm(A,order=-1,err=state)
        call check(error,state%error(),'expm(A, order=-1, err) returns an error')
        if (allocated(error)) return

        B = 1.0_qp
        E = expm(B,err=state)
        call check(error,state%error(),'expm(B, err) of a non-square matrix returns an error')
        if (allocated(error)) return

    end subroutine test_expm_err_form_w

    !> expm without err is pure: it can be called from do concurrent and from a pure procedure
    subroutine test_expm_pure_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        complex(qp) :: A(2,2),E(2,2,3)
        integer :: k

        A = reshape([0.0_qp,1.0_qp,-1.0_qp,0.0_qp], [2,2])

        ! exp(k*A) is a rotation by k radians
        do concurrent(k=1:3)
            E(:,:,k) = expm(k*A)
        end do
        do k = 1,3
            call check(error,maxval(abs(E(:,:,k) - rotation_w(real(k,qp)))) < tol, &
                       'expm in do concurrent')
            if (allocated(error)) return
        end do

        ! exp(2A) = exp(A)**2
        call check(error,maxval(abs(squared_exp_w(A) - E(:,:,2))) < tol,'expm in a pure function')
        if (allocated(error)) return

    end subroutine test_expm_pure_w

    pure function squared_exp_w(A) result(E2)
        complex(qp),intent(in) :: A(:,:)
        complex(qp) :: E2(size(A,1),size(A,2)),E(size(A,1),size(A,2))
        integer :: i,j

        E = expm(A)
        do concurrent(i=1:size(A,1),j=1:size(A,2))
            E2(i,j) = sum(E(i,:)*E(:,j))
        end do

    end function squared_exp_w

    pure function rotation_w(t) result(R)
        real(qp),intent(in) :: t
        complex(qp) :: R(2,2)

        R = reshape([cos(t),sin(t),-sin(t),cos(t)], [2,2])

    end function rotation_w

    !> exp(0) is the identity, and the exponential of an empty matrix is empty
    subroutine test_expm_zero_and_empty_w(error)
        type(error_type),allocatable,intent(out) :: error

        real(qp),parameter :: tol = 1000*epsilon(1.0_qp)
        complex(qp) :: Z(3,3),E(3,3),A0(0,0),E0(0,0)
        complex(qp),allocatable :: F(:,:)
        type(la_state) :: state

        Z = 0.0_qp
        E = expm(Z)
        call check(error,maxval(abs(E - eye(3,mold=1.0_qp))) < tol,'expm of the zero matrix is the identity')
        if (allocated(error)) return

        F = expm(A0)
        call check(error,size(F,1) == 0 .and. size(F,2) == 0,'expm of an empty matrix is empty')
        if (allocated(error)) return

        call matrix_exp(A0,E0,err=state)
        call check(error,state%ok(),'matrix_exp of an empty matrix: '//state%print())
        if (allocated(error)) return

    end subroutine test_expm_zero_and_empty_w
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

end module test_linalg_expm

program test_expm
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_expm,only:test_expm_computation
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_expm",test_expm_computation) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_expm
