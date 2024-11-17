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

        call test_schur_api_s(error)
        if (error) return
        call test_schur_s(error)
        call test_schur_api_d(error)
        if (error) return
        call test_schur_d(error)
        call test_schur_api_q(error)
        if (error) return
        call test_schur_q(error)
        call test_schur_api_c(error)
        if (error) return
        call test_schur_api_z(error)
        if (error) return
        call test_schur_api_w(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Schur tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_schur

    !> QR factorization of a random matrix
    subroutine test_schur_api_s(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp) :: eigs(n)
        real(sp),dimension(n,n) :: a,t,z
        real(sp),allocatable :: storage(:)
        type(linalg_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) return
        
    end subroutine test_schur_api_s
    
    subroutine test_schur_s(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        
        real(sp),dimension(n,n) :: a,t,z
        type(linalg_state) :: err
        
        a = transpose(reshape([[2.65896708,1.42440458,-1.92933439], &
                               [0.,-0.32948354,-0.49063704], &
                               [0.,1.31178921,-0.32948354]], [n,n]))
        
        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return
        
        print "(3(1x,g0.12))",transpose(t)
        print "(3(1x,g0.12))",transpose(z)
        
    end subroutine test_schur_s
    
    subroutine test_schur_api_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp) :: eigs(n)
        real(dp),dimension(n,n) :: a,t,z
        real(dp),allocatable :: storage(:)
        type(linalg_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) return
        
    end subroutine test_schur_api_d
    
    subroutine test_schur_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        
        real(dp),dimension(n,n) :: a,t,z
        type(linalg_state) :: err
        
        a = transpose(reshape([[2.65896708,1.42440458,-1.92933439], &
                               [0.,-0.32948354,-0.49063704], &
                               [0.,1.31178921,-0.32948354]], [n,n]))
        
        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return
        
        print "(3(1x,g0.12))",transpose(t)
        print "(3(1x,g0.12))",transpose(z)
        
    end subroutine test_schur_d
    
    subroutine test_schur_api_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp) :: eigs(n)
        real(qp),dimension(n,n) :: a,t,z
        real(qp),allocatable :: storage(:)
        type(linalg_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) return
        
    end subroutine test_schur_api_q
    
    subroutine test_schur_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        
        real(qp),dimension(n,n) :: a,t,z
        type(linalg_state) :: err
        
        a = transpose(reshape([[2.65896708,1.42440458,-1.92933439], &
                               [0.,-0.32948354,-0.49063704], &
                               [0.,1.31178921,-0.32948354]], [n,n]))
        
        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return
        
        print "(3(1x,g0.12))",transpose(t)
        print "(3(1x,g0.12))",transpose(z)
        
    end subroutine test_schur_q
    
    subroutine test_schur_api_c(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp) :: eigs(n)
        complex(sp),dimension(n,n) :: a,t,z
        complex(sp),allocatable :: storage(:)
        real(sp) :: rea(n,n),ima(n,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) return
        
    end subroutine test_schur_api_c
    
    subroutine test_schur_api_z(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp) :: eigs(n)
        complex(dp),dimension(n,n) :: a,t,z
        complex(dp),allocatable :: storage(:)
        real(dp) :: rea(n,n),ima(n,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) return
        
    end subroutine test_schur_api_z
    
    subroutine test_schur_api_w(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp) :: eigs(n)
        complex(qp),dimension(n,n) :: a,t,z
        complex(qp),allocatable :: storage(:)
        real(qp) :: rea(n,n),ima(n,n)
        type(linalg_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) return
        
    end subroutine test_schur_api_w
    
!import numpy as np
!from scipy.linalg import schur, eigvals
!A = np.array([[0, 2, 2], [0, 1, 2], [1, 0, 1]])
!T, Z = schur(A)
!T
!array([
!Z
!array([[0.72711591, -0.60156188, 0.33079564],
!       [0.52839428, 0.79801892, 0.28976765],
!       [0.43829436, 0.03590414, -0.89811411]])

end module test_linalg_schur

