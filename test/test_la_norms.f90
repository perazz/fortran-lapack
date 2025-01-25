
! Test matrix norms
module test_linalg_norms
    use linear_algebra

    implicit none(type,external)

    contains

    !> Vector norm test: array interfaces
    subroutine test_norms(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_norm_s_1d(error)
        if (error) return
        call test_norm_s_2d(error)
        if (error) return
        call test_norm_s_3d(error)
        if (error) return
        call test_norm_s_4d(error)
        if (error) return
        call test_norm_s_5d(error)
        if (error) return
        call test_norm_s_6d(error)
        if (error) return
        call test_norm2_s_2d(error)
        if (error) return
        call test_norm_dimmed_s_2d(error)
        if (error) return
        call test_norm2_s_3d(error)
        if (error) return
        call test_norm_dimmed_s_3d(error)
        if (error) return
        call test_norm2_s_4d(error)
        if (error) return
        call test_norm_dimmed_s_4d(error)
        if (error) return
        call test_norm2_s_5d(error)
        if (error) return
        call test_norm_dimmed_s_5d(error)
        if (error) return
        call test_norm2_s_6d(error)
        if (error) return
        call test_norm_dimmed_s_6d(error)
        if (error) return
        call test_norm_d_1d(error)
        if (error) return
        call test_norm_d_2d(error)
        if (error) return
        call test_norm_d_3d(error)
        if (error) return
        call test_norm_d_4d(error)
        if (error) return
        call test_norm_d_5d(error)
        if (error) return
        call test_norm_d_6d(error)
        if (error) return
        call test_norm2_d_2d(error)
        if (error) return
        call test_norm_dimmed_d_2d(error)
        if (error) return
        call test_norm2_d_3d(error)
        if (error) return
        call test_norm_dimmed_d_3d(error)
        if (error) return
        call test_norm2_d_4d(error)
        if (error) return
        call test_norm_dimmed_d_4d(error)
        if (error) return
        call test_norm2_d_5d(error)
        if (error) return
        call test_norm_dimmed_d_5d(error)
        if (error) return
        call test_norm2_d_6d(error)
        if (error) return
        call test_norm_dimmed_d_6d(error)
        if (error) return
        call test_norm_q_1d(error)
        if (error) return
        call test_norm_q_2d(error)
        if (error) return
        call test_norm_q_3d(error)
        if (error) return
        call test_norm_q_4d(error)
        if (error) return
        call test_norm_q_5d(error)
        if (error) return
        call test_norm_q_6d(error)
        if (error) return
        call test_norm2_q_2d(error)
        if (error) return
        call test_norm_dimmed_q_2d(error)
        if (error) return
        call test_norm2_q_3d(error)
        if (error) return
        call test_norm_dimmed_q_3d(error)
        if (error) return
        call test_norm2_q_4d(error)
        if (error) return
        call test_norm_dimmed_q_4d(error)
        if (error) return
        call test_norm2_q_5d(error)
        if (error) return
        call test_norm_dimmed_q_5d(error)
        if (error) return
        call test_norm2_q_6d(error)
        if (error) return
        call test_norm_dimmed_q_6d(error)
        if (error) return
        call test_norm_c_1d(error)
        if (error) return
        call test_norm_c_2d(error)
        if (error) return
        call test_norm_c_3d(error)
        if (error) return
        call test_norm_c_4d(error)
        if (error) return
        call test_norm_c_5d(error)
        if (error) return
        call test_norm_c_6d(error)
        if (error) return
        call test_norm_dimmed_c_2d(error)
        if (error) return
        call test_norm_dimmed_c_3d(error)
        if (error) return
        call test_norm_dimmed_c_4d(error)
        if (error) return
        call test_norm_dimmed_c_5d(error)
        if (error) return
        call test_norm_dimmed_c_6d(error)
        if (error) return
        call test_norm_z_1d(error)
        if (error) return
        call test_norm_z_2d(error)
        if (error) return
        call test_norm_z_3d(error)
        if (error) return
        call test_norm_z_4d(error)
        if (error) return
        call test_norm_z_5d(error)
        if (error) return
        call test_norm_z_6d(error)
        if (error) return
        call test_norm_dimmed_z_2d(error)
        if (error) return
        call test_norm_dimmed_z_3d(error)
        if (error) return
        call test_norm_dimmed_z_4d(error)
        if (error) return
        call test_norm_dimmed_z_5d(error)
        if (error) return
        call test_norm_dimmed_z_6d(error)
        if (error) return
        call test_norm_w_1d(error)
        if (error) return
        call test_norm_w_2d(error)
        if (error) return
        call test_norm_w_3d(error)
        if (error) return
        call test_norm_w_4d(error)
        if (error) return
        call test_norm_w_5d(error)
        if (error) return
        call test_norm_w_6d(error)
        if (error) return
        call test_norm_dimmed_w_2d(error)
        if (error) return
        call test_norm_dimmed_w_3d(error)
        if (error) return
        call test_norm_dimmed_w_4d(error)
        if (error) return
        call test_norm_dimmed_w_5d(error)
        if (error) return
        call test_norm_dimmed_w_6d(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Vector norm tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_norms
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_1d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**1
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:)
        
        allocate (a(n),b(2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_s_1d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**2
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_s_2d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**3
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_s_3d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**4
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_s_4d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**5
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_s_5d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**6
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_s_6d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals
    subroutine test_norm2_s_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_s_2d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_s_2d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(sp),allocatable :: a(:),b(:,:)
        real(sp),allocatable :: bnrm(:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_s_2d
    
    subroutine test_norm2_s_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_s_3d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_s_3d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(sp),allocatable :: a(:),b(:,:,:)
        real(sp),allocatable :: bnrm(:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_s_3d
    
    subroutine test_norm2_s_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_s_4d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_s_4d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(sp),allocatable :: a(:),b(:,:,:,:)
        real(sp),allocatable :: bnrm(:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_s_4d
    
    subroutine test_norm2_s_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_s_5d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_s_5d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(sp),allocatable :: a(:),b(:,:,:,:,:)
        real(sp),allocatable :: bnrm(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_s_5d
    
    subroutine test_norm2_s_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_s_6d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_s_6d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:)
        real(sp),allocatable :: bnrm(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_s_6d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_1d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**1
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:)
        
        allocate (a(n),b(2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_d_1d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**2
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_d_2d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**3
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_d_3d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**4
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_d_4d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**5
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_d_5d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**6
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_d_6d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals
    subroutine test_norm2_d_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_d_2d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_d_2d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(dp),allocatable :: a(:),b(:,:)
        real(dp),allocatable :: bnrm(:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_d_2d
    
    subroutine test_norm2_d_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_d_3d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_d_3d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(dp),allocatable :: a(:),b(:,:,:)
        real(dp),allocatable :: bnrm(:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_d_3d
    
    subroutine test_norm2_d_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_d_4d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_d_4d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(dp),allocatable :: a(:),b(:,:,:,:)
        real(dp),allocatable :: bnrm(:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_d_4d
    
    subroutine test_norm2_d_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_d_5d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_d_5d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(dp),allocatable :: a(:),b(:,:,:,:,:)
        real(dp),allocatable :: bnrm(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_d_5d
    
    subroutine test_norm2_d_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_d_6d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_d_6d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:)
        real(dp),allocatable :: bnrm(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_d_6d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_1d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**1
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:)
        
        allocate (a(n),b(2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_q_1d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**2
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_q_2d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**3
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_q_3d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**4
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_q_4d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**5
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_q_5d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**6
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_q_6d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals
    subroutine test_norm2_q_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_q_2d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_q_2d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(qp),allocatable :: a(:),b(:,:)
        real(qp),allocatable :: bnrm(:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_q_2d
    
    subroutine test_norm2_q_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_q_3d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_q_3d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(qp),allocatable :: a(:),b(:,:,:)
        real(qp),allocatable :: bnrm(:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_q_3d
    
    subroutine test_norm2_q_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_q_4d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_q_4d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(qp),allocatable :: a(:),b(:,:,:,:)
        real(qp),allocatable :: bnrm(:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_q_4d
    
    subroutine test_norm2_q_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_q_5d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_q_5d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(qp),allocatable :: a(:),b(:,:,:,:,:)
        real(qp),allocatable :: bnrm(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_q_5d
    
    subroutine test_norm2_q_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) return
            
        end do
        
    end subroutine test_norm2_q_6d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_q_6d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:)
        real(qp),allocatable :: bnrm(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_q_6d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_1d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**1
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:)
        
        allocate (a(n),b(2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_c_1d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**2
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_c_2d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**3
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_c_3d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**4
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_c_4d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**5
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_c_5d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**6
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_sp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_sp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_sp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_c_6d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_c_2d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(sp),allocatable :: a(:),b(:,:)
        real(sp),allocatable :: bnrm(:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_c_2d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_c_3d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(sp),allocatable :: a(:),b(:,:,:)
        real(sp),allocatable :: bnrm(:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_c_3d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_c_4d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(sp),allocatable :: a(:),b(:,:,:,:)
        real(sp),allocatable :: bnrm(:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_c_4d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_c_5d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(sp),allocatable :: a(:),b(:,:,:,:,:)
        real(sp),allocatable :: bnrm(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_c_5d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_c_6d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:)
        real(sp),allocatable :: bnrm(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_c_6d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_1d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**1
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:)
        
        allocate (a(n),b(2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_z_1d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**2
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_z_2d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**3
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_z_3d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**4
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_z_4d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**5
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_z_5d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**6
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_dp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_dp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_dp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_z_6d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_z_2d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(dp),allocatable :: a(:),b(:,:)
        real(dp),allocatable :: bnrm(:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_z_2d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_z_3d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(dp),allocatable :: a(:),b(:,:,:)
        real(dp),allocatable :: bnrm(:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_z_3d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_z_4d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(dp),allocatable :: a(:),b(:,:,:,:)
        real(dp),allocatable :: bnrm(:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_z_4d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_z_5d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(dp),allocatable :: a(:),b(:,:,:,:,:)
        real(dp),allocatable :: bnrm(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_z_5d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_z_6d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:)
        real(dp),allocatable :: bnrm(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_z_6d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_1d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**1
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:)
        
        allocate (a(n),b(2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_w_1d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_2d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**2
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_w_2d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_3d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**3
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_w_3d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_4d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**4
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_w_4d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_5d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**5
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_w_5d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_6d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**6
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        do order = 1,10
           error = .not. abs(norm(a,order) - norm(b,order)) < tol*max(1.0_qp,norm(a,order))
           if (error) return
        end do
        
        ! Infinity norms
        error = .not. abs(norm(a,'inf') - norm(b,'inf')) < tol*max(1.0_qp,norm(a,'inf'))
        if (error) return

        ! Infinity norms
        error = .not. abs(norm(a,'-inf') - norm(b,'-inf')) < tol*max(1.0_qp,norm(a,'-inf'))
        if (error) return
        
    end subroutine test_norm_w_6d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_w_2d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 2
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(qp),allocatable :: a(:),b(:,:)
        real(qp),allocatable :: bnrm(:)
        
        allocate (a(n),b(2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_w_2d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_w_3d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 3
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(qp),allocatable :: a(:),b(:,:,:)
        real(qp),allocatable :: bnrm(:,:)
        
        allocate (a(n),b(2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_w_3d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_w_4d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 4
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(qp),allocatable :: a(:),b(:,:,:,:)
        real(qp),allocatable :: bnrm(:,:,:)
        
        allocate (a(n),b(2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_w_4d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_w_5d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 5
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(qp),allocatable :: a(:),b(:,:,:,:,:)
        real(qp),allocatable :: bnrm(:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_w_5d
    
    ! Test norm along a dimension and compare it against individually evaluated norms
    subroutine test_norm_dimmed_w_6d(error)
        logical,intent(out) :: error
       
        integer(ilp) :: j,dim,order
        integer(ilp),parameter :: ndim = 6
        integer(ilp),parameter :: n = 2_ilp**ndim
        integer(ilp),parameter :: dims(*) = [(dim,dim=1,ndim)]
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        integer(ilp) :: coords(ndim)
        real :: x(ndim)
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:)
        real(qp),allocatable :: bnrm(:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        do order = 1,5
        
           do dim = 1,ndim
            
               bnrm = norm(b,order,dim)
               
               ! Assert size
               error = .not. all(shape(bnrm) == pack(shape(b),dims /= dim))
               if (error) print *, 'INVALID OUTPUT SHAPE, order=',order,' dim=',dim
               
           end do
            
        end do
        
    end subroutine test_norm_dimmed_w_6d
    
end module test_linalg_norms

