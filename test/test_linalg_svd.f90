! Test singular value decomposition
module test_linalg_svd
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> SVD tests
    subroutine test_svd(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_svd_s(error)
        if (error) return
        call test_svd_d(error)
        if (error) return
        call test_svd_q(error)
        if (error) return

        call test_complex_svd_c(error)
        if (error) return
        call test_complex_svd_z(error)
        if (error) return
        call test_complex_svd_w(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('SVD tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_svd

    !> Real matrix svd
    subroutine test_svd_s(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp),parameter :: third = 1.0_sp/3.0_sp
        real(sp),parameter :: twothd = 2*third
        real(sp),parameter :: rsqrt2 = 1.0_sp/sqrt(2.0_sp)
        real(sp),parameter :: rsqrt18 = 1.0_sp/sqrt(18.0_sp)

        real(sp),parameter :: A_mat(2,3) = reshape([real(sp) :: 3,2,2,3,2,-2], [2,3])
        real(sp),parameter :: s_sol(2) = [real(sp) :: 5,3]
        real(sp),parameter :: u_sol(2,2) = reshape(rsqrt2*[1,1,1,-1], [2,2])
        real(sp),parameter :: vt_sol(3,3) = reshape([rsqrt2,rsqrt18,twothd, &
                                                      rsqrt2,-rsqrt18,-twothd, &
                                                      0.0_sp,4*rsqrt18,-third], [3,3])

        !> Local variables
        type(linalg_state) :: state
        real(sp) :: A(2,3),s(2),u(2,2),vt(3,3)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> Function interface
        s = svdvals(A,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> [S, U]. Singular vectors could be all flipped
        call svd(A,s,u,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol))
        if (error) return

        !> [S, U]. Overwrite A matrix
        call svd(A,s,u,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol))
        if (error) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol)) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [S, V^T]. Do not overwrite A matrix
        A = A_mat
        call svd(A,s,vt=vt,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(vt + vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [S, V^T]. Overwrite A matrix
        call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [U, S, V^T].
        A = A_mat
        call svd(A,s,u,vt,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol)) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
        A = A_mat
        u = 0
        vt = 0
        call svd(A,s,u,vt,full_matrices=.false.,err=state)
        error = state%error() &
           .or. .not. all(abs(s - s_sol) <= tol) &
           .or. .not. (all(abs(u(:,:2) - u_sol(:,:2)) <= tol) .or. all(abs(u(:,:2) + u_sol(:,:2)) <= tol)) &
           .or. .not. (all(abs(vt(:2,:) - vt_sol(:2,:)) <= tol) .or. all(abs(vt(:2,:) + vt_sol(:2,:)) <= tol))

        if (error) return

    end subroutine test_svd_s

    subroutine test_svd_d(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp),parameter :: third = 1.0_dp/3.0_dp
        real(dp),parameter :: twothd = 2*third
        real(dp),parameter :: rsqrt2 = 1.0_dp/sqrt(2.0_dp)
        real(dp),parameter :: rsqrt18 = 1.0_dp/sqrt(18.0_dp)

        real(dp),parameter :: A_mat(2,3) = reshape([real(dp) :: 3,2,2,3,2,-2], [2,3])
        real(dp),parameter :: s_sol(2) = [real(dp) :: 5,3]
        real(dp),parameter :: u_sol(2,2) = reshape(rsqrt2*[1,1,1,-1], [2,2])
        real(dp),parameter :: vt_sol(3,3) = reshape([rsqrt2,rsqrt18,twothd, &
                                                      rsqrt2,-rsqrt18,-twothd, &
                                                      0.0_dp,4*rsqrt18,-third], [3,3])

        !> Local variables
        type(linalg_state) :: state
        real(dp) :: A(2,3),s(2),u(2,2),vt(3,3)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> Function interface
        s = svdvals(A,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> [S, U]. Singular vectors could be all flipped
        call svd(A,s,u,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol))
        if (error) return

        !> [S, U]. Overwrite A matrix
        call svd(A,s,u,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol))
        if (error) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol)) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [S, V^T]. Do not overwrite A matrix
        A = A_mat
        call svd(A,s,vt=vt,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(vt + vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [S, V^T]. Overwrite A matrix
        call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [U, S, V^T].
        A = A_mat
        call svd(A,s,u,vt,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol)) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
        A = A_mat
        u = 0
        vt = 0
        call svd(A,s,u,vt,full_matrices=.false.,err=state)
        error = state%error() &
           .or. .not. all(abs(s - s_sol) <= tol) &
           .or. .not. (all(abs(u(:,:2) - u_sol(:,:2)) <= tol) .or. all(abs(u(:,:2) + u_sol(:,:2)) <= tol)) &
           .or. .not. (all(abs(vt(:2,:) - vt_sol(:2,:)) <= tol) .or. all(abs(vt(:2,:) + vt_sol(:2,:)) <= tol))

        if (error) return

    end subroutine test_svd_d

    subroutine test_svd_q(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp),parameter :: third = 1.0_qp/3.0_qp
        real(qp),parameter :: twothd = 2*third
        real(qp),parameter :: rsqrt2 = 1.0_qp/sqrt(2.0_qp)
        real(qp),parameter :: rsqrt18 = 1.0_qp/sqrt(18.0_qp)

        real(qp),parameter :: A_mat(2,3) = reshape([real(qp) :: 3,2,2,3,2,-2], [2,3])
        real(qp),parameter :: s_sol(2) = [real(qp) :: 5,3]
        real(qp),parameter :: u_sol(2,2) = reshape(rsqrt2*[1,1,1,-1], [2,2])
        real(qp),parameter :: vt_sol(3,3) = reshape([rsqrt2,rsqrt18,twothd, &
                                                      rsqrt2,-rsqrt18,-twothd, &
                                                      0.0_qp,4*rsqrt18,-third], [3,3])

        !> Local variables
        type(linalg_state) :: state
        real(qp) :: A(2,3),s(2),u(2,2),vt(3,3)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> Function interface
        s = svdvals(A,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> [S, U]. Singular vectors could be all flipped
        call svd(A,s,u,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol))
        if (error) return

        !> [S, U]. Overwrite A matrix
        call svd(A,s,u,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol))
        if (error) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol)) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [S, V^T]. Do not overwrite A matrix
        A = A_mat
        call svd(A,s,vt=vt,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(vt + vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [S, V^T]. Overwrite A matrix
        call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [U, S, V^T].
        A = A_mat
        call svd(A,s,u,vt,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. (all(abs(u - u_sol) <= tol) .or. all(abs(u + u_sol) <= tol)) .or. &
                .not. (all(abs(vt - vt_sol) <= tol) .or. all(abs(vt + vt_sol) <= tol))
        if (error) return

        !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
        A = A_mat
        u = 0
        vt = 0
        call svd(A,s,u,vt,full_matrices=.false.,err=state)
        error = state%error() &
           .or. .not. all(abs(s - s_sol) <= tol) &
           .or. .not. (all(abs(u(:,:2) - u_sol(:,:2)) <= tol) .or. all(abs(u(:,:2) + u_sol(:,:2)) <= tol)) &
           .or. .not. (all(abs(vt(:2,:) - vt_sol(:2,:)) <= tol) .or. all(abs(vt(:2,:) + vt_sol(:2,:)) <= tol))

        if (error) return

    end subroutine test_svd_q

    !> Test complex svd
    subroutine test_complex_svd_c(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp),parameter :: one = 1.0_sp
        real(sp),parameter :: zero = 0.0_sp
        real(sp),parameter :: sqrt2 = sqrt(2.0_sp)
        real(sp),parameter :: rsqrt2 = one/sqrt2
        complex(sp),parameter :: cone = (1.0_sp,0.0_sp)
        complex(sp),parameter :: cimg = (0.0_sp,1.0_sp)
        complex(sp),parameter :: czero = (0.0_sp,0.0_sp)

        real(sp),parameter :: s_sol(2) = [sqrt2,sqrt2]
        complex(sp),parameter :: A_mat(2,2) = reshape([cone,cimg,cimg,cone], [2,2])
        complex(sp),parameter :: u_sol(2,2) = reshape(rsqrt2*[cone,cimg,cimg,cone], [2,2])
        complex(sp),parameter :: vt_sol(2,2) = reshape([cone,czero,czero,cone], [2,2])

        !> Local variables
        type(linalg_state) :: state
        complex(sp) :: A(2,2),u(2,2),vt(2,2)
        real(sp) :: s(2)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> Function interface
        s = svdvals(A,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. all(abs(matmul(u,matmul(diag(s),vt)) - A_mat) <= tol)
        if (error) return

    end subroutine test_complex_svd_c

    subroutine test_complex_svd_z(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp),parameter :: one = 1.0_dp
        real(dp),parameter :: zero = 0.0_dp
        real(dp),parameter :: sqrt2 = sqrt(2.0_dp)
        real(dp),parameter :: rsqrt2 = one/sqrt2
        complex(dp),parameter :: cone = (1.0_dp,0.0_dp)
        complex(dp),parameter :: cimg = (0.0_dp,1.0_dp)
        complex(dp),parameter :: czero = (0.0_dp,0.0_dp)

        real(dp),parameter :: s_sol(2) = [sqrt2,sqrt2]
        complex(dp),parameter :: A_mat(2,2) = reshape([cone,cimg,cimg,cone], [2,2])
        complex(dp),parameter :: u_sol(2,2) = reshape(rsqrt2*[cone,cimg,cimg,cone], [2,2])
        complex(dp),parameter :: vt_sol(2,2) = reshape([cone,czero,czero,cone], [2,2])

        !> Local variables
        type(linalg_state) :: state
        complex(dp) :: A(2,2),u(2,2),vt(2,2)
        real(dp) :: s(2)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> Function interface
        s = svdvals(A,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. all(abs(matmul(u,matmul(diag(s),vt)) - A_mat) <= tol)
        if (error) return

    end subroutine test_complex_svd_z

    subroutine test_complex_svd_w(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp),parameter :: one = 1.0_qp
        real(qp),parameter :: zero = 0.0_qp
        real(qp),parameter :: sqrt2 = sqrt(2.0_qp)
        real(qp),parameter :: rsqrt2 = one/sqrt2
        complex(qp),parameter :: cone = (1.0_qp,0.0_qp)
        complex(qp),parameter :: cimg = (0.0_qp,1.0_qp)
        complex(qp),parameter :: czero = (0.0_qp,0.0_qp)

        real(qp),parameter :: s_sol(2) = [sqrt2,sqrt2]
        complex(qp),parameter :: A_mat(2,2) = reshape([cone,cimg,cimg,cone], [2,2])
        complex(qp),parameter :: u_sol(2,2) = reshape(rsqrt2*[cone,cimg,cimg,cone], [2,2])
        complex(qp),parameter :: vt_sol(2,2) = reshape([cone,czero,czero,cone], [2,2])

        !> Local variables
        type(linalg_state) :: state
        complex(qp) :: A(2,2),u(2,2),vt(2,2)
        real(qp) :: s(2)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> Function interface
        s = svdvals(A,err=state)
        error = state%error() .or. .not. all(abs(s - s_sol) <= tol)
        if (error) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        error = state%error() .or. &
                .not. all(abs(s - s_sol) <= tol) .or. &
                .not. all(abs(matmul(u,matmul(diag(s),vt)) - A_mat) <= tol)
        if (error) return

    end subroutine test_complex_svd_w

end module test_linalg_svd

