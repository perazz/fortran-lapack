! Test QR factorization
module test_linalg_qr
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> SVD tests
    subroutine test_qr(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_qr_random_s(error)
        if (error) return
        call test_qr_random_d(error)
        if (error) return
        call test_qr_random_q(error)
        if (error) return
        call test_qr_random_c(error)
        if (error) return
        call test_qr_random_z(error)
        if (error) return
        call test_qr_random_w(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('QR tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_qr

    !> QR factorization of a random matrix
    subroutine test_qr_random_s(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(sp) :: rea(m,n),ima(m,n)
        integer(ilp) :: lwork
        real(sp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = 0.0_sp
        q = 1.0_sp/q
        r = 0.0_sp
        r = 1.0_sp/r
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 3) overwrite A
        call qr(a,qred,rred,overwrite_a=.true.,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 4) External storage option
        a = rea
        call qr_space(a,lwork)
        allocate (work(lwork))
        call qr(a,q,r,storage=work,err=err)
    
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
    end subroutine test_qr_random_s

    subroutine test_qr_random_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(dp) :: rea(m,n),ima(m,n)
        integer(ilp) :: lwork
        real(dp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = 0.0_dp
        q = 1.0_dp/q
        r = 0.0_dp
        r = 1.0_dp/r
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 3) overwrite A
        call qr(a,qred,rred,overwrite_a=.true.,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 4) External storage option
        a = rea
        call qr_space(a,lwork)
        allocate (work(lwork))
        call qr(a,q,r,storage=work,err=err)
    
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
    end subroutine test_qr_random_d

    subroutine test_qr_random_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(qp) :: rea(m,n),ima(m,n)
        integer(ilp) :: lwork
        real(qp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = 0.0_qp
        q = 1.0_qp/q
        r = 0.0_qp
        r = 1.0_qp/r
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 3) overwrite A
        call qr(a,qred,rred,overwrite_a=.true.,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 4) External storage option
        a = rea
        call qr_space(a,lwork)
        allocate (work(lwork))
        call qr(a,q,r,storage=work,err=err)
    
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
    end subroutine test_qr_random_q

    subroutine test_qr_random_c(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(sp) :: rea(m,n),ima(m,n)
        integer(ilp) :: lwork
        complex(sp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = 0.0_sp
        q = 1.0_sp/q
        r = 0.0_sp
        r = 1.0_sp/r
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 3) overwrite A
        call qr(a,qred,rred,overwrite_a=.true.,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 4) External storage option
        a = cmplx(rea,ima,kind=sp)
        call qr_space(a,lwork)
        allocate (work(lwork))
        call qr(a,q,r,storage=work,err=err)
    
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
    end subroutine test_qr_random_c

    subroutine test_qr_random_z(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(dp) :: rea(m,n),ima(m,n)
        integer(ilp) :: lwork
        complex(dp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = 0.0_dp
        q = 1.0_dp/q
        r = 0.0_dp
        r = 1.0_dp/r
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 3) overwrite A
        call qr(a,qred,rred,overwrite_a=.true.,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 4) External storage option
        a = cmplx(rea,ima,kind=dp)
        call qr_space(a,lwork)
        allocate (work(lwork))
        call qr(a,q,r,storage=work,err=err)
    
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
    end subroutine test_qr_random_z

    subroutine test_qr_random_w(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(qp) :: rea(m,n),ima(m,n)
        integer(ilp) :: lwork
        complex(qp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        
        ! 1) QR factorization with full matrices. Input NaNs to be sure Q and R are OK on return
        q = 0.0_qp
        q = 1.0_qp/q
        r = 0.0_qp
        r = 1.0_qp/r
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 3) overwrite A
        call qr(a,qred,rred,overwrite_a=.true.,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < tol)
        if (error) return
        
        ! 4) External storage option
        a = cmplx(rea,ima,kind=qp)
        call qr_space(a,lwork)
        allocate (work(lwork))
        call qr(a,q,r,storage=work,err=err)
    
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < tol)
        if (error) return
        
    end subroutine test_qr_random_w

end module test_linalg_qr

