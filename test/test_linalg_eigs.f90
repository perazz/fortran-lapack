! Test eigendecomposition
module test_linalg_eig
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> SVD tests
    subroutine test_eig(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_eig_diagonal_s(error)
        if (error) return
        call test_eig_diagonal_d(error)
        if (error) return
        call test_eig_diagonal_q(error)
        if (error) return
        
        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Eigenproblem tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_eig

    !> Simple diagonal matrix eigenvalues
    subroutine test_eig_diagonal_s(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: zero = 0.0_sp
        real(sp),parameter :: tol = sqrt(epsilon(zero))

        !> Local variables
        type(linalg_state) :: state
        real(sp) :: A(3,3),B(2,2)
        complex(sp) :: lambda(3),lambdab(2),Bvec(2,2)

        !> Matrix with real eigenvalues
        A = reshape([1,0,0, &
                     0,2,0, &
                     0,0,3], [3,3])
        
        call eig(A,lambda,err=state)
        error = state%error() .or. &
                .not. all(aimag(lambda) == zero .and. real(lambda,kind=sp) == [1,2,3])
        if (error) return
        
        !> Matrix with complex eigenvalues
        B = transpose(reshape([1,-1, &
                               1,1], [2,2]))
        
        call eig(B,lambda,left=Bvec,err=state)
        error = state%error()
        if (error) return
        
!eigenvalues, eigenvectors = LA.eig(np.array([[1, -1], [1, 1]]))
!eigenvalues
!array([1.+1.j, 1.-1.j])
!eigenvectors
!array([[0.70710678+0.j        , 0.70710678-0.j        ],
!       [0.        -0.70710678j, 0.        +0.70710678j]])
!

    end subroutine test_eig_diagonal_s

    subroutine test_eig_diagonal_d(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: zero = 0.0_dp
        real(dp),parameter :: tol = sqrt(epsilon(zero))

        !> Local variables
        type(linalg_state) :: state
        real(dp) :: A(3,3),B(2,2)
        complex(dp) :: lambda(3),lambdab(2),Bvec(2,2)

        !> Matrix with real eigenvalues
        A = reshape([1,0,0, &
                     0,2,0, &
                     0,0,3], [3,3])
        
        call eig(A,lambda,err=state)
        error = state%error() .or. &
                .not. all(aimag(lambda) == zero .and. real(lambda,kind=dp) == [1,2,3])
        if (error) return
        
        !> Matrix with complex eigenvalues
        B = transpose(reshape([1,-1, &
                               1,1], [2,2]))
        
        call eig(B,lambda,left=Bvec,err=state)
        error = state%error()
        if (error) return
        
!eigenvalues, eigenvectors = LA.eig(np.array([[1, -1], [1, 1]]))
!eigenvalues
!array([1.+1.j, 1.-1.j])
!eigenvectors
!array([[0.70710678+0.j        , 0.70710678-0.j        ],
!       [0.        -0.70710678j, 0.        +0.70710678j]])
!

    end subroutine test_eig_diagonal_d

    subroutine test_eig_diagonal_q(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: zero = 0.0_qp
        real(qp),parameter :: tol = sqrt(epsilon(zero))

        !> Local variables
        type(linalg_state) :: state
        real(qp) :: A(3,3),B(2,2)
        complex(qp) :: lambda(3),lambdab(2),Bvec(2,2)

        !> Matrix with real eigenvalues
        A = reshape([1,0,0, &
                     0,2,0, &
                     0,0,3], [3,3])
        
        call eig(A,lambda,err=state)
        error = state%error() .or. &
                .not. all(aimag(lambda) == zero .and. real(lambda,kind=qp) == [1,2,3])
        if (error) return
        
        !> Matrix with complex eigenvalues
        B = transpose(reshape([1,-1, &
                               1,1], [2,2]))
        
        call eig(B,lambda,left=Bvec,err=state)
        error = state%error()
        if (error) return
        
!eigenvalues, eigenvectors = LA.eig(np.array([[1, -1], [1, 1]]))
!eigenvalues
!array([1.+1.j, 1.-1.j])
!eigenvectors
!array([[0.70710678+0.j        , 0.70710678-0.j        ],
!       [0.        -0.70710678j, 0.        +0.70710678j]])
!

    end subroutine test_eig_diagonal_q

end module test_linalg_eig

