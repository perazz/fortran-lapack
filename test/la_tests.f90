program la_tests
    use test_la_aux
    use test_la_eye
    use test_la_solve
    use test_la_inverse
    use test_la_least_squares
    use test_la_determinant
    use test_la_svd
    use test_la_eigs
    use test_la_qr
    use test_la_norms
    use test_la_schur
    use test_la_pinv
    implicit none(type, external)

    integer :: i,seed_size
    integer, allocatable :: seed(:)
    logical :: error

    !> Fixed seed: the suite must draw the same random matrices on every run
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = [(1000 + 7*i,i=1,seed_size)]
    call random_seed(put=seed)

    call test_formats(error)
    if (error) error stop 'test_formats'

    call test_solve(error)
    if (error) error stop 'test_solve'

    call test_inverse_matrix(error)
    if (error) error stop 'test_inverse_matrix'

    call test_least_squares(error)
    if (error) error stop 'test_least_squares'

    call test_matrix_determinant(error)
    if (error) error stop 'test_determinant'

    call test_eye(error)
    if (error) error stop 'test_eye'

    call test_svd(error)
    if (error) error stop 'test_svd'
    
    call test_eig(error)
    if (error) error stop 'test_eig'

    call test_norms(error)
    if (error) error stop 'test_norms'
    
    call test_schur(error)
    if (error) error stop 'test_schur'

    call test_pseudoinverse_matrix(error)
    if (error) error stop 'test_pseudoinverse'
    
    !> All tests passed
    stop 0

end program la_tests
