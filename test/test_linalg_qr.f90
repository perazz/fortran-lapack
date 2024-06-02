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

        call test_qr_s(error)
        if (error) return
        call test_qr_d(error)
        if (error) return
        call test_qr_q(error)
        if (error) return
        call test_qr_c(error)
        if (error) return
        call test_qr_z(error)
        if (error) return
        call test_qr_w(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('QR tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_qr

    !> QR factorization of a random matrix
    subroutine test_qr_s(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(sp) :: a(m,n),q(m,k),r(k,n)
        real(sp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_sp)))
        if (error) return
        
    end subroutine test_qr_s

    subroutine test_qr_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(dp) :: a(m,n),q(m,k),r(k,n)
        real(dp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_dp)))
        if (error) return
        
    end subroutine test_qr_d

    subroutine test_qr_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        real(qp) :: a(m,n),q(m,k),r(k,n)
        real(qp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        a = rea
        
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_qp)))
        if (error) return
        
    end subroutine test_qr_q

    subroutine test_qr_c(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        complex(sp) :: a(m,n),q(m,k),r(k,n)
        real(sp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_sp)))
        if (error) return
        
    end subroutine test_qr_c

    subroutine test_qr_z(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        complex(dp) :: a(m,n),q(m,k),r(k,n)
        real(dp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_dp)))
        if (error) return
        
    end subroutine test_qr_z

    subroutine test_qr_w(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: m = 15_ilp
        integer(ilp),parameter :: n = 4_ilp
        integer(ilp),parameter :: k = min(m,n)
        complex(qp) :: a(m,n),q(m,k),r(k,n)
        real(qp) :: rea(m,n),ima(m,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        
        call qr(a,q,r,err=err)
        
        ! Check return code
        error = err%error()
        if (error) return
        
        ! Check solution
        error = .not. all(abs(a - matmul(q,r)) < sqrt(epsilon(0.0_qp)))
        if (error) return
        
    end subroutine test_qr_w

end module test_linalg_qr

