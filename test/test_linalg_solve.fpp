#:set REAL_KINDS    = ["sp", "dp", "qp"]
#:set REAL_INITIALS = ["s","d","q"]
#:set REAL_TYPES    = ["real({})".format(k) for k in REAL_KINDS]
#:set REAL_KINDS_TYPES = list(zip(REAL_KINDS, REAL_TYPES, REAL_INITIALS))

module test_linalg_solve
    use stdlib_linalg_constants
    use stdlib_linalg_state
    use stdlib_linalg_solve

    implicit none (type,external)

    contains

    !> Solve real linear system
    subroutine test_solve(error)
        logical, intent(out) :: error

        #:for rk,rt,ri in REAL_KINDS_TYPES
        call test_${ri}$solve(error)
        if (error) return
        #: endfor

    end subroutine test_solve

    !> Simple linear system
    #:for rk,rt,ri in REAL_KINDS_TYPES
    subroutine test_${ri}$solve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        ${rt}$ :: A(3,3) = transpose(reshape([${rt}$ :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        ${rt}$ :: b  (3) = [${rt}$ :: 1, 4, -1]
        ${rt}$ :: res(3) = [${rt}$ :: -2, -2, 3]
        ${rt}$ :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_${rk}$)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_${rk}$)
        print *, 'state = ',state%print()


    end subroutine test_${ri}$solve

    #:endfor

end module test_linalg_solve


