! Test singular value decomposition
module test_linalg_svd
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:diag,svd,svdvals
    use la_state_type,only:la_state
    implicit none(type,external)
    
    public :: test_svd

    contains

    !> Solve several SVD problems
    subroutine test_svd(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)
        
        allocate (tests(0))

        call add_test(tests,new_unittest("test_svd_s",test_svd_s))
        call add_test(tests,new_unittest("test_svd_d",test_svd_d))
        call add_test(tests,new_unittest("test_svd_q",test_svd_q))

        call add_test(tests,new_unittest("test_complex_svd_c",test_complex_svd_c))
        call add_test(tests,new_unittest("test_complex_svd_z",test_complex_svd_z))
        call add_test(tests,new_unittest("test_complex_svd_w",test_complex_svd_w))

    end subroutine test_svd

    !> Real matrix svd
    subroutine test_svd_s(error)
        type(error_type),allocatable,intent(out) :: error

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
        character(:),allocatable :: test
        type(la_state) :: state
        real(sp) :: A(2,3),s(2),u(2,2),vt(3,3)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        
        test = 'subroutine version'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        
        !> Function interface
        s = svdvals(A,err=state)
        
        test = 'function interface'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return

        !> [S, U]. Singular vectors could be all flipped
        call svd(A,s,u,err=state)
        
        test = 'subroutine with singular vectors'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return

        !> [S, U]. Overwrite A matrix
        call svd(A,s,u,overwrite_a=.true.,err=state)
        
        test = 'subroutine, overwrite_a'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        
        test = '[S, U, V^T]'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return
        
        !> [S, V^T]. Do not overwrite A matrix
        A = A_mat
        call svd(A,s,vt=vt,err=state)

        test = '[S, V^T], overwrite_a=.false.'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return

        !> [S, V^T]. Overwrite A matrix
        call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
        
        test = '[S, V^T], overwrite_a=.true.'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return
        
        !> [U, S, V^T].
        A = A_mat
        call svd(A,s,u,vt,err=state)
        
        test = '[U, S, V^T]'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return

        !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
        A = A_mat
        u = 0
        vt = 0
        call svd(A,s,u,vt,full_matrices=.false.,err=state)
        
        test = '[U, S, V^T], partial storage'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(abs(u(:,:2)) - abs(u_sol(:,:2))) <= tol),test//': U(:,:2)')
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt(:2,:)) - abs(vt_sol(:2,:))) <= tol),test//': V^T(:2,:)')
        if (allocated(error)) return

    end subroutine test_svd_s

    subroutine test_svd_d(error)
        type(error_type),allocatable,intent(out) :: error

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
        character(:),allocatable :: test
        type(la_state) :: state
        real(dp) :: A(2,3),s(2),u(2,2),vt(3,3)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        
        test = 'subroutine version'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        
        !> Function interface
        s = svdvals(A,err=state)
        
        test = 'function interface'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return

        !> [S, U]. Singular vectors could be all flipped
        call svd(A,s,u,err=state)
        
        test = 'subroutine with singular vectors'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return

        !> [S, U]. Overwrite A matrix
        call svd(A,s,u,overwrite_a=.true.,err=state)
        
        test = 'subroutine, overwrite_a'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        
        test = '[S, U, V^T]'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return
        
        !> [S, V^T]. Do not overwrite A matrix
        A = A_mat
        call svd(A,s,vt=vt,err=state)

        test = '[S, V^T], overwrite_a=.false.'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return

        !> [S, V^T]. Overwrite A matrix
        call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
        
        test = '[S, V^T], overwrite_a=.true.'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return
        
        !> [U, S, V^T].
        A = A_mat
        call svd(A,s,u,vt,err=state)
        
        test = '[U, S, V^T]'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return

        !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
        A = A_mat
        u = 0
        vt = 0
        call svd(A,s,u,vt,full_matrices=.false.,err=state)
        
        test = '[U, S, V^T], partial storage'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(abs(u(:,:2)) - abs(u_sol(:,:2))) <= tol),test//': U(:,:2)')
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt(:2,:)) - abs(vt_sol(:2,:))) <= tol),test//': V^T(:2,:)')
        if (allocated(error)) return

    end subroutine test_svd_d

    subroutine test_svd_q(error)
        type(error_type),allocatable,intent(out) :: error

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
        character(:),allocatable :: test
        type(la_state) :: state
        real(qp) :: A(2,3),s(2),u(2,2),vt(3,3)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        
        test = 'subroutine version'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        
        !> Function interface
        s = svdvals(A,err=state)
        
        test = 'function interface'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return

        !> [S, U]. Singular vectors could be all flipped
        call svd(A,s,u,err=state)
        
        test = 'subroutine with singular vectors'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return

        !> [S, U]. Overwrite A matrix
        call svd(A,s,u,overwrite_a=.true.,err=state)
        
        test = 'subroutine, overwrite_a'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        
        test = '[S, U, V^T]'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return
        
        !> [S, V^T]. Do not overwrite A matrix
        A = A_mat
        call svd(A,s,vt=vt,err=state)

        test = '[S, V^T], overwrite_a=.false.'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return

        !> [S, V^T]. Overwrite A matrix
        call svd(A,s,vt=vt,overwrite_a=.true.,err=state)
        
        test = '[S, V^T], overwrite_a=.true.'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return
        
        !> [U, S, V^T].
        A = A_mat
        call svd(A,s,u,vt,err=state)
        
        test = '[U, S, V^T]'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(abs(u) - abs(u_sol)) <= tol),test//': U')
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt) - abs(vt_sol)) <= tol),test//': V^T')
        if (allocated(error)) return

        !> [U, S, V^T]. Partial storage -> compare until k=2 columns of U rows of V^T
        A = A_mat
        u = 0
        vt = 0
        call svd(A,s,u,vt,full_matrices=.false.,err=state)
        
        test = '[U, S, V^T], partial storage'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(abs(u(:,:2)) - abs(u_sol(:,:2))) <= tol),test//': U(:,:2)')
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(abs(vt(:2,:)) - abs(vt_sol(:2,:))) <= tol),test//': V^T(:2,:)')
        if (allocated(error)) return

    end subroutine test_svd_q

    !> Test complex svd
    subroutine test_complex_svd_c(error)
        type(error_type),allocatable,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
        real(sp),parameter :: one = 1.0_sp
        real(sp),parameter :: zero = 0.0_sp
        real(sp),parameter :: sqrt2 = sqrt(2.0_sp)
        real(sp),parameter :: rsqrt2 = one/sqrt2
        complex(sp),parameter :: csqrt2 = (rsqrt2,zero)
        complex(sp),parameter :: isqrt2 = (zero,rsqrt2)
        complex(sp),parameter :: cone = (1.0_sp,0.0_sp)
        complex(sp),parameter :: cimg = (0.0_sp,1.0_sp)
        complex(sp),parameter :: czero = (0.0_sp,0.0_sp)

        real(sp),parameter :: s_sol(2) = [sqrt2,sqrt2]
        complex(sp),parameter :: A_mat(2,2) = reshape([cone,cimg,cimg,cone], [2,2])
        complex(sp),parameter :: u_sol(2,2) = reshape([csqrt2,isqrt2,isqrt2,csqrt2], [2,2])
        complex(sp),parameter :: vt_sol(2,2) = reshape([cone,czero,czero,cone], [2,2])

        !> Local variables
        character(:),allocatable :: test
        type(la_state) :: state
        complex(sp) :: A(2,2),u(2,2),vt(2,2)
        real(sp) :: s(2)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        
        test = '[S], complex subroutine'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        
        !> Function interface
        s = svdvals(A,err=state)

        test = 'svdvals, complex function'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        
        test = '[S, U, V^T], complex'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(matmul(u,matmul(diag(s),vt)) - A_mat) <= tol),test//': U*S*V^T')
        if (allocated(error)) return

    end subroutine test_complex_svd_c

    subroutine test_complex_svd_z(error)
        type(error_type),allocatable,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
        real(dp),parameter :: one = 1.0_dp
        real(dp),parameter :: zero = 0.0_dp
        real(dp),parameter :: sqrt2 = sqrt(2.0_dp)
        real(dp),parameter :: rsqrt2 = one/sqrt2
        complex(dp),parameter :: csqrt2 = (rsqrt2,zero)
        complex(dp),parameter :: isqrt2 = (zero,rsqrt2)
        complex(dp),parameter :: cone = (1.0_dp,0.0_dp)
        complex(dp),parameter :: cimg = (0.0_dp,1.0_dp)
        complex(dp),parameter :: czero = (0.0_dp,0.0_dp)

        real(dp),parameter :: s_sol(2) = [sqrt2,sqrt2]
        complex(dp),parameter :: A_mat(2,2) = reshape([cone,cimg,cimg,cone], [2,2])
        complex(dp),parameter :: u_sol(2,2) = reshape([csqrt2,isqrt2,isqrt2,csqrt2], [2,2])
        complex(dp),parameter :: vt_sol(2,2) = reshape([cone,czero,czero,cone], [2,2])

        !> Local variables
        character(:),allocatable :: test
        type(la_state) :: state
        complex(dp) :: A(2,2),u(2,2),vt(2,2)
        real(dp) :: s(2)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        
        test = '[S], complex subroutine'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        
        !> Function interface
        s = svdvals(A,err=state)

        test = 'svdvals, complex function'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        
        test = '[S, U, V^T], complex'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(matmul(u,matmul(diag(s),vt)) - A_mat) <= tol),test//': U*S*V^T')
        if (allocated(error)) return

    end subroutine test_complex_svd_z

    subroutine test_complex_svd_w(error)
        type(error_type),allocatable,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
        real(qp),parameter :: one = 1.0_qp
        real(qp),parameter :: zero = 0.0_qp
        real(qp),parameter :: sqrt2 = sqrt(2.0_qp)
        real(qp),parameter :: rsqrt2 = one/sqrt2
        complex(qp),parameter :: csqrt2 = (rsqrt2,zero)
        complex(qp),parameter :: isqrt2 = (zero,rsqrt2)
        complex(qp),parameter :: cone = (1.0_qp,0.0_qp)
        complex(qp),parameter :: cimg = (0.0_qp,1.0_qp)
        complex(qp),parameter :: czero = (0.0_qp,0.0_qp)

        real(qp),parameter :: s_sol(2) = [sqrt2,sqrt2]
        complex(qp),parameter :: A_mat(2,2) = reshape([cone,cimg,cimg,cone], [2,2])
        complex(qp),parameter :: u_sol(2,2) = reshape([csqrt2,isqrt2,isqrt2,csqrt2], [2,2])
        complex(qp),parameter :: vt_sol(2,2) = reshape([cone,czero,czero,cone], [2,2])

        !> Local variables
        character(:),allocatable :: test
        type(la_state) :: state
        complex(qp) :: A(2,2),u(2,2),vt(2,2)
        real(qp) :: s(2)

        !> Initialize matrix
        A = A_mat

        !> Simple subroutine version
        call svd(A,s,err=state)
        
        test = '[S], complex subroutine'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        
        !> Function interface
        s = svdvals(A,err=state)

        test = 'svdvals, complex function'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return

        !> [S, U, V^T]
        A = A_mat
        call svd(A,s,u,vt,overwrite_a=.true.,err=state)
        
        test = '[S, U, V^T], complex'
        call check(error,state%ok(),test//': '//state%print())
        if (allocated(error)) return
        call check(error,all(abs(s - s_sol) <= tol),test//': S')
        if (allocated(error)) return
        call check(error,all(abs(matmul(u,matmul(diag(s),vt)) - A_mat) <= tol),test//': U*S*V^T')
        if (allocated(error)) return

    end subroutine test_complex_svd_w

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type),allocatable,intent(inout) :: tests(:)
        type(unittest_type),intent(in) :: new_test
        
        integer :: n
        type(unittest_type),allocatable :: new_tests(:)
        
        if (allocated(tests)) then
            n = size(tests)
        else
            n = 0
        end if
        
        allocate (new_tests(n + 1))
        if (n > 0) new_tests(1:n) = tests(1:n)
                 new_tests(1 + n) = new_test
        call move_alloc(from=new_tests,to=tests)
    
    end subroutine add_test

end module test_linalg_svd

program test_svd_driver
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_svd,only:test_svd
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_svd",test_svd) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_svd_driver
