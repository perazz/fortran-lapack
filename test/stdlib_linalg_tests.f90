program la_linalg_tests
    use test_linalg_aux
    use test_linalg_eye
    use test_linalg_solve
    use test_linalg_inverse
    use test_linalg_least_squares
    use test_linalg_determinant
    use test_linalg_svd
    use test_linalg_eig
    use test_linalg_qr
    use test_linalg_norms
    use test_linalg_schur
    use test_linalg_pseudoinverse
    implicit none(type, external)

    logical :: error

!    call test_formats(error)
!    if (error) error stop 'test_formats'
!
!    call test_solve(error)
!    if (error) error stop 'test_solve'
!
!    call test_inverse_matrix(error)
!    if (error) error stop 'test_inverse_matrix'
!
!    call test_least_squares(error)
!    if (error) error stop 'test_least_squares'
!
!    call test_matrix_determinant(error)
!    if (error) error stop 'test_determinant'
!
!    call test_eye(error)
!    if (error) error stop 'test_eye'
!
!    call test_svd(error)
!    if (error) error stop 'test_svd'
!    
!    call test_eig(error)
!    if (error) error stop 'test_eig'

!    call test_norms(error)
!    if (error) error stop 'test_norms'
!    call test_schur(error)
!    if (error) error stop 'test_schur'

    call test_pseudoinverse_matrix(error)
    if (error) error stop 'test_pseudoinverse'
    
    !> All tests passed
    stop 0

end program la_linalg_tests
