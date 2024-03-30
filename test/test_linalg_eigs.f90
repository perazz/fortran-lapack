! Test eigendecomposition
module test_linalg_eigs
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> SVD tests
    subroutine test_eig(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_eig_s(error)
        if (error) return
        call test_eig_d(error)
        if (error) return
        call test_eig_q(error)
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

    end subroutine test_eig

    !> Real matrix svd
    subroutine test_eig_s(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))

        !> Local variables
        type(linalg_state) :: state
        real(sp) :: A(3,3)
        complex(sp) :: lambda

        !> Initialize matrix
        A = diag([1,2,3])
        
        call eig(A,lambda,err=state)
        error = state%error()

        if (error) return

    end subroutine test_eig_s

    subroutine test_eig_d(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))

        !> Local variables
        type(linalg_state) :: state
        real(dp) :: A(3,3)
        complex(dp) :: lambda

        !> Initialize matrix
        A = diag([1,2,3])
        
        call eig(A,lambda,err=state)
        error = state%error()

        if (error) return

    end subroutine test_eig_d

    subroutine test_eig_q(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))

        !> Local variables
        type(linalg_state) :: state
        real(qp) :: A(3,3)
        complex(qp) :: lambda

        !> Initialize matrix
        A = diag([1,2,3])
        
        call eig(A,lambda,err=state)
        error = state%error()

        if (error) return

    end subroutine test_eig_q

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

end module test_linalg_eigs

