program stdlib_linalg_tests
    use test_linalg_aux
    use test_linalg_solve
    use test_linalg_inverse
    implicit none(type, external)

    logical :: error


    call test_formats(error)
    if (error) error stop 'test_formats'

    call test_solve(error)
    if (error) error stop 'test_solve'

    call test_inverse_matrix(error)
    if (error) error stop 'test_inverse_matrix'

    !> All tests passed
    stop 0


end program stdlib_linalg_tests
