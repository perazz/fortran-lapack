! Test eigendecomposition
module test_linalg_eig
    use la_interface

    implicit none(type,external)

    contains

    !> SVD tests
    subroutine test_eig(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_eig_real_s(error)
        if (error) return
        call test_eig_real_d(error)
        if (error) return
        call test_eig_real_q(error)
        if (error) return
        
        call test_eigh_real_s(error)
        if (error) return
        call test_eigh_real_d(error)
        if (error) return
        call test_eigh_real_q(error)
        if (error) return
                
        call test_eig_complex_c(error)
        if (error) return
        call test_eig_complex_z(error)
        if (error) return
        call test_eig_complex_w(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Eigenproblem tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_eig

    !> Simple real matrix eigenvalues
    subroutine test_eig_real_s(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: zero = 0.0_sp
        real(sp),parameter :: two = 2.0_sp
        real(sp),parameter :: sqrt2o2 = sqrt(two)*0.5_sp
        real(sp),parameter :: tol = sqrt(epsilon(zero))

        !> Local variables
        type(linalg_state) :: state
        real(sp) :: A(3,3),B(2,2)
        complex(sp) :: lambda(3),Bvec(2,2),Bres(2,2)

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
                               
        !> Expected right eigenvectors
        Bres(1,1:2) = sqrt2o2
        Bres(2,1) = cmplx(zero,-sqrt2o2,kind=sp)
        Bres(2,2) = cmplx(zero,+sqrt2o2,kind=sp)
        
        call eig(B,lambda,right=Bvec,err=state)
        error = state%error() .or. any(abs(Bres - Bvec) > tol)
        
        print *, bvec(1,:)
        print *, bvec(2,:)
        print *, bres
        
        if (error) return
        
    end subroutine test_eig_real_s

    ! Symmetric matrix eigenvalues
    subroutine test_eigh_real_s(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: zero = 0.0_sp
        real(sp),parameter :: tol = sqrt(epsilon(zero))
        real(sp),parameter :: A(4,4) = reshape([6,3,1,5, &
                                               3,0,5,1, &
                                               1,5,6,2, &
                                               5,1,2,2], [4,4])
        
        !> Local variables
        real(sp) :: Amat(4,4),lambda(4),vect(4,4),Av(4,4),lv(4,4)
        type(linalg_state) :: state
        
        Amat = A
        
        call eigh(Amat,lambda,vect,err=state)
        
        Av = matmul(A,vect)
        lv = matmul(vect,diag(lambda))
        
        error = state%error() .or. .not. all(abs(Av - lv) < tol*abs(Av))
        if (error) return
        
        !> Test functional versions
        lambda = eigvalsh(Amat)
   
        lambda = eigvalsh(Amat,err=state)
        error = state%error()
        if (error) return
        
        !> Test functional versions
        Amat = A
        lambda = eigvalsh(Amat,upper_a=.false.,err=state)
        error = state%error()
        if (error) return
        
    end subroutine test_eigh_real_s

    subroutine test_eig_real_d(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: zero = 0.0_dp
        real(dp),parameter :: two = 2.0_dp
        real(dp),parameter :: sqrt2o2 = sqrt(two)*0.5_dp
        real(dp),parameter :: tol = sqrt(epsilon(zero))

        !> Local variables
        type(linalg_state) :: state
        real(dp) :: A(3,3),B(2,2)
        complex(dp) :: lambda(3),Bvec(2,2),Bres(2,2)

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
                               
        !> Expected right eigenvectors
        Bres(1,1:2) = sqrt2o2
        Bres(2,1) = cmplx(zero,-sqrt2o2,kind=dp)
        Bres(2,2) = cmplx(zero,+sqrt2o2,kind=dp)
        
        call eig(B,lambda,right=Bvec,err=state)
        error = state%error() .or. any(abs(Bres - Bvec) > tol)
        
        print *, bvec(1,:)
        print *, bvec(2,:)
        print *, bres
        
        if (error) return
        
    end subroutine test_eig_real_d

    ! Symmetric matrix eigenvalues
    subroutine test_eigh_real_d(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: zero = 0.0_dp
        real(dp),parameter :: tol = sqrt(epsilon(zero))
        real(dp),parameter :: A(4,4) = reshape([6,3,1,5, &
                                               3,0,5,1, &
                                               1,5,6,2, &
                                               5,1,2,2], [4,4])
        
        !> Local variables
        real(dp) :: Amat(4,4),lambda(4),vect(4,4),Av(4,4),lv(4,4)
        type(linalg_state) :: state
        
        Amat = A
        
        call eigh(Amat,lambda,vect,err=state)
        
        Av = matmul(A,vect)
        lv = matmul(vect,diag(lambda))
        
        error = state%error() .or. .not. all(abs(Av - lv) < tol*abs(Av))
        if (error) return
        
        !> Test functional versions
        lambda = eigvalsh(Amat)
   
        lambda = eigvalsh(Amat,err=state)
        error = state%error()
        if (error) return
        
        !> Test functional versions
        Amat = A
        lambda = eigvalsh(Amat,upper_a=.false.,err=state)
        error = state%error()
        if (error) return
        
    end subroutine test_eigh_real_d

    subroutine test_eig_real_q(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: zero = 0.0_qp
        real(qp),parameter :: two = 2.0_qp
        real(qp),parameter :: sqrt2o2 = sqrt(two)*0.5_qp
        real(qp),parameter :: tol = sqrt(epsilon(zero))

        !> Local variables
        type(linalg_state) :: state
        real(qp) :: A(3,3),B(2,2)
        complex(qp) :: lambda(3),Bvec(2,2),Bres(2,2)

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
                               
        !> Expected right eigenvectors
        Bres(1,1:2) = sqrt2o2
        Bres(2,1) = cmplx(zero,-sqrt2o2,kind=qp)
        Bres(2,2) = cmplx(zero,+sqrt2o2,kind=qp)
        
        call eig(B,lambda,right=Bvec,err=state)
        error = state%error() .or. any(abs(Bres - Bvec) > tol)
        
        print *, bvec(1,:)
        print *, bvec(2,:)
        print *, bres
        
        if (error) return
        
    end subroutine test_eig_real_q

    ! Symmetric matrix eigenvalues
    subroutine test_eigh_real_q(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: zero = 0.0_qp
        real(qp),parameter :: tol = sqrt(epsilon(zero))
        real(qp),parameter :: A(4,4) = reshape([6,3,1,5, &
                                               3,0,5,1, &
                                               1,5,6,2, &
                                               5,1,2,2], [4,4])
        
        !> Local variables
        real(qp) :: Amat(4,4),lambda(4),vect(4,4),Av(4,4),lv(4,4)
        type(linalg_state) :: state
        
        Amat = A
        
        call eigh(Amat,lambda,vect,err=state)
        
        Av = matmul(A,vect)
        lv = matmul(vect,diag(lambda))
        
        error = state%error() .or. .not. all(abs(Av - lv) < tol*abs(Av))
        if (error) return
        
        !> Test functional versions
        lambda = eigvalsh(Amat)
   
        lambda = eigvalsh(Amat,err=state)
        error = state%error()
        if (error) return
        
        !> Test functional versions
        Amat = A
        lambda = eigvalsh(Amat,upper_a=.false.,err=state)
        error = state%error()
        if (error) return
        
    end subroutine test_eigh_real_q

    !> Simple complex matrix eigenvalues
    subroutine test_eig_complex_c(error)
        logical,intent(out) :: error

        !> Reference solution
        real(sp),parameter :: zero = 0.0_sp
        real(sp),parameter :: two = 2.0_sp
        real(sp),parameter :: sqrt2o2 = sqrt(two)*0.5_sp
        real(sp),parameter :: tol = sqrt(epsilon(zero))
        complex(sp),parameter :: cone = (1.0_sp,0.0_sp)
        complex(sp),parameter :: cimg = (0.0_sp,1.0_sp)
        complex(sp),parameter :: czero = (0.0_sp,0.0_sp)

        !> Local vaciables
        type(linalg_state) :: state
        complex(sp) :: A(2,2),lambda(2),Avec(2,2),Ares(2,2),lres(2)

        !> Matcix with real eigenvalues
        A = transpose(reshape([cone,cimg, &
                               -cimg,cone], [2,2]))
                
        call eig(A,lambda,right=Avec,err=state)
        
        !> Expected eigenvalues and eigenvectors
        lres(1) = two
        lres(2) = zero
        
        !> Eigenvectors may vary: do not use for error
        Ares(1,1) = cmplx(zero,sqrt2o2,kind=sp)
        Ares(1,2) = cmplx(sqrt2o2,zero,kind=sp)
        Ares(2,1) = cmplx(sqrt2o2,zero,kind=sp)
        Ares(2,2) = cmplx(zero,sqrt2o2,kind=sp)
        
        error = state%error() .or. any(abs(lambda - lres) > tol)
        if (error) return
        
    end subroutine test_eig_complex_c

    subroutine test_eig_complex_z(error)
        logical,intent(out) :: error

        !> Reference solution
        real(dp),parameter :: zero = 0.0_dp
        real(dp),parameter :: two = 2.0_dp
        real(dp),parameter :: sqrt2o2 = sqrt(two)*0.5_dp
        real(dp),parameter :: tol = sqrt(epsilon(zero))
        complex(dp),parameter :: cone = (1.0_dp,0.0_dp)
        complex(dp),parameter :: cimg = (0.0_dp,1.0_dp)
        complex(dp),parameter :: czero = (0.0_dp,0.0_dp)

        !> Local vaciables
        type(linalg_state) :: state
        complex(dp) :: A(2,2),lambda(2),Avec(2,2),Ares(2,2),lres(2)

        !> Matcix with real eigenvalues
        A = transpose(reshape([cone,cimg, &
                               -cimg,cone], [2,2]))
                
        call eig(A,lambda,right=Avec,err=state)
        
        !> Expected eigenvalues and eigenvectors
        lres(1) = two
        lres(2) = zero
        
        !> Eigenvectors may vary: do not use for error
        Ares(1,1) = cmplx(zero,sqrt2o2,kind=dp)
        Ares(1,2) = cmplx(sqrt2o2,zero,kind=dp)
        Ares(2,1) = cmplx(sqrt2o2,zero,kind=dp)
        Ares(2,2) = cmplx(zero,sqrt2o2,kind=dp)
        
        error = state%error() .or. any(abs(lambda - lres) > tol)
        if (error) return
        
    end subroutine test_eig_complex_z

    subroutine test_eig_complex_w(error)
        logical,intent(out) :: error

        !> Reference solution
        real(qp),parameter :: zero = 0.0_qp
        real(qp),parameter :: two = 2.0_qp
        real(qp),parameter :: sqrt2o2 = sqrt(two)*0.5_qp
        real(qp),parameter :: tol = sqrt(epsilon(zero))
        complex(qp),parameter :: cone = (1.0_qp,0.0_qp)
        complex(qp),parameter :: cimg = (0.0_qp,1.0_qp)
        complex(qp),parameter :: czero = (0.0_qp,0.0_qp)

        !> Local vaciables
        type(linalg_state) :: state
        complex(qp) :: A(2,2),lambda(2),Avec(2,2),Ares(2,2),lres(2)

        !> Matcix with real eigenvalues
        A = transpose(reshape([cone,cimg, &
                               -cimg,cone], [2,2]))
                
        call eig(A,lambda,right=Avec,err=state)
        
        !> Expected eigenvalues and eigenvectors
        lres(1) = two
        lres(2) = zero
        
        !> Eigenvectors may vary: do not use for error
        Ares(1,1) = cmplx(zero,sqrt2o2,kind=qp)
        Ares(1,2) = cmplx(sqrt2o2,zero,kind=qp)
        Ares(2,1) = cmplx(sqrt2o2,zero,kind=qp)
        Ares(2,2) = cmplx(zero,sqrt2o2,kind=qp)
        
        error = state%error() .or. any(abs(lambda - lres) > tol)
        if (error) return
        
    end subroutine test_eig_complex_w

end module test_linalg_eig

