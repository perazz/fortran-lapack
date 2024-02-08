
module test_linalg_solve
    use stdlib_linalg_constants
    use stdlib_linalg_state
    use stdlib_linalg_solve

    implicit none (type,external)

    contains

    !> Solve real linear system
    subroutine test_solve(error)
        logical, intent(out) :: error

        call test_ssolve(error)
        if (error) return
        call test_dsolve(error)
        if (error) return
        call test_qsolve(error)
        if (error) return

    end subroutine test_solve

    !> Simple linear system
    subroutine test_ssolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(sp) :: A(3,3) = transpose(reshape([real(sp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(sp) :: b  (3) = [real(sp) :: 1, 4, -1]
        real(sp) :: res(3) = [real(sp) :: -2, -2, 3]
        real(sp) :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_sp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_sp)
        print *, 'state = ',state%print()


    end subroutine test_ssolve

    subroutine test_dsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(dp) :: A(3,3) = transpose(reshape([real(dp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(dp) :: b  (3) = [real(dp) :: 1, 4, -1]
        real(dp) :: res(3) = [real(dp) :: -2, -2, 3]
        real(dp) :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_dp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_dp)
        print *, 'state = ',state%print()


    end subroutine test_dsolve

    subroutine test_qsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(qp) :: A(3,3) = transpose(reshape([real(qp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(qp) :: b  (3) = [real(qp) :: 1, 4, -1]
        real(qp) :: res(3) = [real(qp) :: -2, -2, 3]
        real(qp) :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_qp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_qp)
        print *, 'state = ',state%print()


    end subroutine test_qsolve


end module test_linalg_solve


