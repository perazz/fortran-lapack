! Test Schur form
module test_la_schur
    use linear_algebra

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
#ifdef LA_WITH_XDP
        call test_schur_api_x(error)
        if (error) return
        call test_schur_x(error)
#endif
#ifdef LA_WITH_QP
        call test_schur_api_q(error)
        if (error) return
        call test_schur_q(error)
#endif
        call test_schur_api_c(error)
        if (error) return
        call test_schur_api_z(error)
        if (error) return
#ifdef LA_WITH_XDP
        call test_schur_api_y(error)
        if (error) return
#endif
#ifdef LA_WITH_QP
        call test_schur_api_w(error)
        if (error) return
#endif

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
        type(la_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_s
    
    subroutine test_schur_s(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(sp),parameter :: rtol = 1.0e-4_sp
        real(sp),parameter :: eps = sqrt(epsilon(0.0_sp))
        real(sp),dimension(n,n) :: a,t,z
        type(la_state) :: err

        call random_number(a)

        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        error = .not. (all(abs(a - matmul(matmul(z,t),transpose(z))) <= max(rtol*abs(z),eps)))
        if (error) print *, 'invalid matmul real(sp)'
        if (error) print *, a - matmul(matmul(z,t),transpose(z))
                
    end subroutine test_schur_s
    
    subroutine test_schur_api_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp) :: eigs(n)
        real(dp),dimension(n,n) :: a,t,z
        real(dp),allocatable :: storage(:)
        type(la_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_d
    
    subroutine test_schur_d(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(dp),parameter :: rtol = 1.0e-4_dp
        real(dp),parameter :: eps = sqrt(epsilon(0.0_dp))
        real(dp),dimension(n,n) :: a,t,z
        type(la_state) :: err

        call random_number(a)

        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        error = .not. (all(abs(a - matmul(matmul(z,t),transpose(z))) <= max(rtol*abs(z),eps)))
        if (error) print *, 'invalid matmul real(dp)'
        if (error) print *, a - matmul(matmul(z,t),transpose(z))
                
    end subroutine test_schur_d
    
#ifdef LA_WITH_XDP
    subroutine test_schur_api_x(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(xdp),parameter :: tol = 10*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: eigs(n)
        real(xdp),dimension(n,n) :: a,t,z
        real(xdp),allocatable :: storage(:)
        type(la_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_x
    
    subroutine test_schur_x(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(xdp),parameter :: rtol = 1.0e-4_xdp
        real(xdp),parameter :: eps = sqrt(epsilon(0.0_xdp))
        real(xdp),dimension(n,n) :: a,t,z
        type(la_state) :: err

        call random_number(a)

        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        error = .not. (all(abs(a - matmul(matmul(z,t),transpose(z))) <= max(rtol*abs(z),eps)))
        if (error) print *, 'invalid matmul real(xdp)'
        if (error) print *, a - matmul(matmul(z,t),transpose(z))
                
    end subroutine test_schur_x
    
#endif
    
#ifdef LA_WITH_QP
    subroutine test_schur_api_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp) :: eigs(n)
        real(qp),dimension(n,n) :: a,t,z
        real(qp),allocatable :: storage(:)
        type(la_state) :: err
        
        call random_number(a)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_q
    
    subroutine test_schur_q(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp
        real(qp),parameter :: rtol = 1.0e-4_qp
        real(qp),parameter :: eps = sqrt(epsilon(0.0_qp))
        real(qp),dimension(n,n) :: a,t,z
        type(la_state) :: err

        call random_number(a)

        ! Run schur
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        error = .not. (all(abs(a - matmul(matmul(z,t),transpose(z))) <= max(rtol*abs(z),eps)))
        if (error) print *, 'invalid matmul real(qp)'
        if (error) print *, a - matmul(matmul(z,t),transpose(z))
                
    end subroutine test_schur_q
    
#endif
    
    subroutine test_schur_api_c(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp) :: eigs(n)
        complex(sp),dimension(n,n) :: a,t,z
        complex(sp),allocatable :: storage(:)
        real(sp) :: rea(n,n),ima(n,n)
        type(la_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
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
        type(la_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_z
    
#ifdef LA_WITH_XDP
    subroutine test_schur_api_y(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(xdp),parameter :: tol = 10*sqrt(epsilon(0.0_xdp))
        complex(xdp) :: eigs(n)
        complex(xdp),dimension(n,n) :: a,t,z
        complex(xdp),allocatable :: storage(:)
        real(xdp) :: rea(n,n),ima(n,n)
        type(la_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=xdp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_y
    
#endif
    
#ifdef LA_WITH_QP
    subroutine test_schur_api_w(error)
        logical,intent(out) :: error

        integer(ilp),parameter :: n = 15_ilp
        integer(ilp) :: lwork
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp) :: eigs(n)
        complex(qp),dimension(n,n) :: a,t,z
        complex(qp),allocatable :: storage(:)
        real(qp) :: rea(n,n),ima(n,n)
        type(la_state) :: err
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=qp)
        
        ! Test simple API
        call schur(a,t,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test output transformation matrix
        call schur(a,t,z,err=err)
        error = err%error(); if (error) print err%print(); if (error) return

        ! Test output eigenvalues
        call schur(a,t,eigvals=eigs,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test storage query
        call schur_space(a,lwork,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
        ! Test with user-defined storage
        allocate (storage(lwork))
        call schur(a,t,eigvals=eigs,storage=storage,err=err)
        error = err%error(); if (error) print err%print(); if (error) return
        
    end subroutine test_schur_api_w
    
#endif
    
end module test_la_schur

