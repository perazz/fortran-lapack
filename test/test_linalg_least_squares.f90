module test_linalg_least_squares
    use stdlib_linalg_constants
    use stdlib_linalg_state
    use stdlib_linalg_least_squares

    implicit none (type,external)

    contains

    !> Solve sample least square problems
    subroutine test_least_squares(error)
        logical, intent(out) :: error

        print *, 'Least squares tests...'

        call test_slstsq_one(error)
        if (error) return
        call test_dlstsq_one(error)
        if (error) return
        call test_qlstsq_one(error)
        if (error) return

    end subroutine test_least_squares

    !> Simple polynomial fit
    subroutine test_slstsq_one(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        !> Example scattered data
        real(sp), parameter :: x(*)  = real([1.0, 2.5, 3.5, 4.0, 5.0, 7.0, 8.5], sp)
        real(sp), parameter :: y(*)  = real([0.3, 1.1, 1.5, 2.0, 3.2, 6.6, 8.6], sp)
        real(sp), parameter :: ab(*) = real([0.20925829,  0.12013861], sp)

        real(sp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,err=state)

        print *, 'real(sp): p = ',p

        error = state%error() .or. .not.all(abs(p-ab)<1.0e-6_sp)

    end subroutine test_slstsq_one
    subroutine test_dlstsq_one(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        !> Example scattered data
        real(dp), parameter :: x(*)  = real([1.0, 2.5, 3.5, 4.0, 5.0, 7.0, 8.5], dp)
        real(dp), parameter :: y(*)  = real([0.3, 1.1, 1.5, 2.0, 3.2, 6.6, 8.6], dp)
        real(dp), parameter :: ab(*) = real([0.20925829,  0.12013861], dp)

        real(dp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,err=state)

        print *, 'real(dp): p = ',p

        error = state%error() .or. .not.all(abs(p-ab)<1.0e-6_dp)

    end subroutine test_dlstsq_one
    subroutine test_qlstsq_one(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        !> Example scattered data
        real(qp), parameter :: x(*)  = real([1.0, 2.5, 3.5, 4.0, 5.0, 7.0, 8.5], qp)
        real(qp), parameter :: y(*)  = real([0.3, 1.1, 1.5, 2.0, 3.2, 6.6, 8.6], qp)
        real(qp), parameter :: ab(*) = real([0.20925829,  0.12013861], qp)

        real(qp) :: M(size(x),2),p(2)

        ! Coefficient matrix for polynomial y = a + b*x**2
        M(:,1) = x**0
        M(:,2) = x**2

        ! Find polynomial
        p = lstsq(M,y,err=state)

        print *, 'real(qp): p = ',p

        error = state%error() .or. .not.all(abs(p-ab)<1.0e-6_qp)

    end subroutine test_qlstsq_one

end module test_linalg_least_squares


