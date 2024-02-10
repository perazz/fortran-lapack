program stdlib_linalg_tests
    use test_linalg_aux
    use test_linalg_solve
    implicit none(type, external)

    logical :: error


    call test_formats(error)
    if (error) error stop 'test_formats'

    call test_solve(error)
    if (error) error stop 'test_solve'

    !> All tests passed
    stop 0


end program stdlib_linalg_tests
