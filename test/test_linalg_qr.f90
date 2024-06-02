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
        real(sp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(sp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        ! 1) QR factorization with full matrices
        call qr(a,q,r,err=err)
        print *, 'full, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        print *, 'epsilon=',epsilon(0.0_sp),' within bounds=', &
        count(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_sp))),' out of ',m*n
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_sp)))
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        print *, 'reduced, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        print *, 'epsilon=',epsilon(0.0_sp),' within bounds=', &
        count(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_sp))),' out of ',m*n
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_sp)))
        
        if (error) return
    
    end subroutine test_qr_random_s

    subroutine test_qr_random_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(dp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        ! 1) QR factorization with full matrices
        call qr(a,q,r,err=err)
        print *, 'full, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        print *, 'epsilon=',epsilon(0.0_dp),' within bounds=', &
        count(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_dp))),' out of ',m*n
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_dp)))
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        print *, 'reduced, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        print *, 'epsilon=',epsilon(0.0_dp),' within bounds=', &
        count(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_dp))),' out of ',m*n
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_dp)))
        
        if (error) return
    
    end subroutine test_qr_random_d

    subroutine test_qr_random_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(qp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        ! 1) QR factorization with full matrices
        call qr(a,q,r,err=err)
        print *, 'full, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        print *, 'epsilon=',epsilon(0.0_qp),' within bounds=', &
        count(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_qp))),' out of ',m*n
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_qp)))
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        print *, 'reduced, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        print *, 'epsilon=',epsilon(0.0_qp),' within bounds=', &
        count(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_qp))),' out of ',m*n
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_qp)))
        
        if (error) return
    
    end subroutine test_qr_random_q

    subroutine test_qr_random_c(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        complex(sp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(sp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! 1) QR factorization with full matrices
        call qr(a,q,r,err=err)
        print *, 'full, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        print *, 'epsilon=',epsilon(0.0_sp),' within bounds=', &
        count(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_sp))),' out of ',m*n
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_sp)))
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        print *, 'reduced, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        print *, 'epsilon=',epsilon(0.0_sp),' within bounds=', &
        count(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_sp))),' out of ',m*n
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_sp)))
        
        if (error) return
    
    end subroutine test_qr_random_c

    subroutine test_qr_random_z(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        complex(dp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(dp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! 1) QR factorization with full matrices
        call qr(a,q,r,err=err)
        print *, 'full, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        print *, 'epsilon=',epsilon(0.0_dp),' within bounds=', &
        count(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_dp))),' out of ',m*n
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_dp)))
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        print *, 'reduced, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        print *, 'epsilon=',epsilon(0.0_dp),' within bounds=', &
        count(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_dp))),' out of ',m*n
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_dp)))
        
        if (error) return
    
    end subroutine test_qr_random_z

    subroutine test_qr_random_w(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        complex(qp) :: a(m,n),q(m,m),r(m,n),qred(m,k),rred(k,n)
        real(qp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        
        ! 1) QR factorization with full matrices
        call qr(a,q,r,err=err)
        print *, 'full, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        print *, 'epsilon=',epsilon(0.0_qp),' within bounds=', &
        count(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_qp))),' out of ',m*n
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_qp)))
        if (error) return
        
        ! 2) QR factorization with reduced matrices
        call qr(a,qred,rred,err=err)
        print *, 'reduced, err=',err%print()
        
        ! Check return code
        error = err%error()
        if (error) return
        
        print *, 'epsilon=',epsilon(0.0_qp),' within bounds=', &
        count(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_qp))),' out of ',m*n
        
        ! Check solution
        error = .not. all(abs(a - matmul(qred,rred)) < sqrt(epsilon(0.0_qp)))
        
        if (error) return
    
    end subroutine test_qr_random_w

end module test_linalg_qr

