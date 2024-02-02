program stdlib_linalg_tests
    use test_linalg_aux
    implicit none(type, external)

    logical :: error

    call test_formats(error)
    if (error) error stop 'test_formats'



    !> All tests passed
    stop 0


end program stdlib_linalg_tests
