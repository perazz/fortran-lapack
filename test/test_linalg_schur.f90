! Test Schur form
module test_linalg_schur
    use stdlib_linalg_interface

    implicit none(type,external)

    contains

    !> SVD tests
    subroutine test_schur(error)
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

1       format('Schur tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_schur

    !> QR factorization of a random matrix
    subroutine test_qr_random_s(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp) :: a(n,n),t(n,n)
        integer(ilp) :: lwork
        real(sp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(a)
        
        ! 1) QR factorization with full matrices
        call schur(a,t,err=err)
        
        print *, err%print()
        
        error = err%error()
        
        ! Check solution
        if (error) return
        
    end subroutine test_qr_random_s

    subroutine test_qr_random_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp) :: a(n,n),t(n,n)
        integer(ilp) :: lwork
        real(dp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(a)
        
        ! 1) QR factorization with full matrices
        call schur(a,t,err=err)
        
        print *, err%print()
        
        error = err%error()
        
        ! Check solution
        if (error) return
        
    end subroutine test_qr_random_d

    subroutine test_qr_random_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp) :: a(n,n),t(n,n)
        integer(ilp) :: lwork
        real(qp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(a)
        
        ! 1) QR factorization with full matrices
        call schur(a,t,err=err)
        
        print *, err%print()
        
        error = err%error()
        
        ! Check solution
        if (error) return
        
    end subroutine test_qr_random_q

    subroutine test_qr_random_c(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp) :: a(n,n),t(n,n)
        real(sp) :: rea(n,n),ima(n,n)
        integer(ilp) :: lwork
        complex(sp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! 1) QR factorization with full matrices
        call schur(a,t,err=err)
        
        print *, err%print()
        
        error = err%error()
        
        ! Check solution
        if (error) return
        
    end subroutine test_qr_random_c

    subroutine test_qr_random_z(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp) :: a(n,n),t(n,n)
        real(dp) :: rea(n,n),ima(n,n)
        integer(ilp) :: lwork
        complex(dp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! 1) QR factorization with full matrices
        call schur(a,t,err=err)
        
        print *, err%print()
        
        error = err%error()
        
        ! Check solution
        if (error) return
        
    end subroutine test_qr_random_z

    subroutine test_qr_random_w(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp) :: a(n,n),t(n,n)
        real(qp) :: rea(n,n),ima(n,n)
        integer(ilp) :: lwork
        complex(qp),allocatable :: work(:)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        
        ! 1) QR factorization with full matrices
        call schur(a,t,err=err)
        
        print *, err%print()
        
        error = err%error()
        
        ! Check solution
        if (error) return
        
    end subroutine test_qr_random_w

end module test_linalg_schur

