program stdlib_linalg_tests
    use test_linalg_aux
    use test_linalg_eye
    use test_linalg_solve
    use test_linalg_inverse
    use test_linalg_least_squares
    use test_linalg_determinant
    implicit none(type, external)

    logical :: error
!
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

    call test_svd
    stop 'after test_svd'

    !> All tests passed
    stop 0

    contains


        subroutine test_svd()

            real(dp), parameter :: third   = 1.0_dp/3.0_dp
            real(dp), parameter :: twothd  = 2*third
            real(dp), parameter :: rsqrt2  = 1.0_dp/sqrt(2.0_dp)
            real(dp), parameter :: rsqrt18 = 1.0_dp/sqrt(18.0_dp)

            real(dp) :: A(2,3),s(2),u(2,2),vt(3,3)
            real(dp), parameter :: s_values(2) = [real(dp) :: 5, 3]

            real(dp), parameter ::  u_matrix(2,2) = reshape(rsqrt2*[1,1,1,-1],[2,2])
            real(dp), parameter :: vt_matrix(3,3) = reshape([rsqrt2,rsqrt18,twothd, &
                                                             rsqrt2,-rsqrt18,-twothd,&
                                                             0.0_dp,4*rsqrt18,-third],[3,3])


            A = reshape([real(dp) :: 3,2, 2,3, 2,-2],[2,3])

            call svd(A,s)
            print *, '   s only   '
            print *, 'A=',A
            print *, 's=',s

            call svd(A,s,u)
            print *, '   s,u only   '
            print *, 'A=',A
            print *, 's=',s
            print *, 'u=',u

            call svd(A,s,u,vt=vt)
            print *, '   s,u,vt   '
            print *, 'A=',A
            print *, 's=',s
            print *, 'u=',u
            print *, 'vt=',vt

        end subroutine test_svd


end program stdlib_linalg_tests
