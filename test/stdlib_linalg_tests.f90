program stdlib_linalg_tests
    use test_linalg_aux
    use test_linalg_eye
    use test_linalg_solve
    use test_linalg_inverse
    use test_linalg_least_squares
    use test_linalg_determinant
    use test_linalg_svd
    implicit none(type, external)

    logical :: error

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

    !> All tests passed
    stop 0

!    contains
!
!
!        !> Test real SVD
!        subroutine test_svd(error)
!            logical,intent(out) :: error
!
!            !> Reference solution
!            real(dp), parameter :: tol     = sqrt(epsilon(0.0_dp))
!            real(dp), parameter :: third   = 1.0_dp/3.0_dp
!            real(dp), parameter :: twothd  = 2*third
!            real(dp), parameter :: rsqrt2  = 1.0_dp/sqrt(2.0_dp)
!            real(dp), parameter :: rsqrt18 = 1.0_dp/sqrt(18.0_dp)
!
!            real(dp), parameter ::  A_mat(2,3) = reshape([real(dp) :: 3,2, 2,3, 2,-2],[2,3])
!            real(dp), parameter ::  s_sol(2)   = [real(dp) :: 5, 3]
!            real(dp), parameter ::  u_sol(2,2) = reshape(rsqrt2*[1,1,1,-1],[2,2])
!            real(dp), parameter :: vt_sol(3,3) = reshape([rsqrt2,rsqrt18,twothd, &
!                                                          rsqrt2,-rsqrt18,-twothd,&
!                                                          0.0_dp,4*rsqrt18,-third],[3,3])
!
!            !> Local variables
!            type(linalg_state) :: state
!            real(dp) :: A(2,3),s(2),u(2,2),vt(3,3)
!
!            !> Initialize matrix
!            A = A_mat
!
!            !> Simple subroutine version
!            call svd(A,s,err=state)
!            error = state%error() .or. .not. all(abs(s-s_sol)<=tol)
!            if (error) return
!
!            !> Function interface
!            s = svdvals(A,err=state)
!            error = state%error() .or. .not. all(abs(s-s_sol)<=tol)
!            if (error) return
!
!            !> [S, U]. Singular vectors could be all flipped
!            call svd(A,s,u,err=state)
!            error = state%error() .or. &
!                    .not. all(abs(s-s_sol)<=tol) .or. &
!                    .not.(all(abs(u-u_sol)<=tol) .or. all(abs(u+u_sol)<=tol))
!            if (error) return
!
!            !> [S, U]. Overwrite A matrix
!            call svd(A,s,u,overwrite_a=.true.,err=state)
!            error = state%error() .or. &
!                    .not. all(abs(s-s_sol)<=tol) .or. &
!                    .not.(all(abs(u-u_sol)<=tol) .or. all(abs(u+u_sol)<=tol))
!            if (error) return
!
!            !> [S, U, V^T]
!            A = A_mat
!            call svd(A,s,u,vt,overwrite_a=.true.,err=state)
!            error = state%error() .or. &
!                    .not. all(abs(s-s_sol)<=tol) .or. &
!                    .not.(all(abs(u-u_sol)<=tol) .or. all(abs(u+u_sol)<=tol)) .or. &
!                    .not.(all(abs(vt-vt_sol)<=tol) .or. all(abs(vt+vt_sol)<=tol))
!            if (error) return
!
!            !> [S, V^T]. Do not overwrite A matrix
!            A = A_mat
!            call svd(A,s,vt=vt,err=state)
!            error = state%error() .or. &
!                    .not. all(abs(s-s_sol)<=tol) .or. &
!                    .not.(all(abs(vt+vt_sol)<=tol) .or. all(abs(vt+vt_sol)<=tol))
!            if (error) return
!
!            !> [S, V^T]. Overwrite A matrix
!            call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
!            error = state%error() .or. &
!                    .not. all(abs(s-s_sol)<=tol) .or. &
!                    .not.(all(abs(vt-vt_sol)<=tol) .or. all(abs(vt+vt_sol)<=tol))
!            if (error) return
!
!            !> [U, S, V^T].
!            A = A_mat
!            call svd(A,s,u,vt,err=state)
!            error = state%error() .or. &
!                    .not. all(abs(s-s_sol)<=tol) .or. &
!                    .not.(all(abs(u-u_sol)<=tol) .or. all(abs(u+u_sol)<=tol)) .or. &
!                    .not.(all(abs(vt-vt_sol)<=tol) .or. all(abs(vt+vt_sol)<=tol))
!            if (error) return
!
!            !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
!            A  = A_mat
!            u  = 0
!            vt = 0
!            call svd(A,s,u,vt,full_matrices=.false.,err=state)
!            error = state%error() &
!               .or. .not. all(abs(s-s_sol)<=tol) &
!               .or. .not.(all(abs( u(:,:2)- u_sol(:,:2))<=tol) .or. all(abs( u(:,:2)+ u_sol(:,:2))<=tol)) &
!               .or. .not.(all(abs(vt(:2,:)-vt_sol(:2,:))<=tol) .or. all(abs(vt(:2,:)+vt_sol(:2,:))<=tol))
!
!            if (error) return
!
!        end subroutine test_svd


end program stdlib_linalg_tests
