
! Test matrix norms
module test_linalg_norms
    use stdlib_linalg_interface

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
        call test_norm_s_7d(error)
        if (error) return
        call test_norm_s_8d(error)
        if (error) return
        call test_norm_s_9d(error)
        if (error) return
        call test_norm_s_10d(error)
        if (error) return
        call test_norm_s_11d(error)
        if (error) return
        call test_norm_s_12d(error)
        if (error) return
        call test_norm_s_13d(error)
        if (error) return
        call test_norm_s_14d(error)
        if (error) return
        call test_norm2_s_2d(error)
        if (error) return
        call test_norm2_s_3d(error)
        if (error) return
        call test_norm2_s_4d(error)
        if (error) return
        call test_norm2_s_5d(error)
        if (error) return
        call test_norm2_s_6d(error)
        if (error) return
        call test_norm2_s_7d(error)
        if (error) return
        call test_norm2_s_8d(error)
        if (error) return
        call test_norm2_s_9d(error)
        if (error) return
        call test_norm2_s_10d(error)
        if (error) return
        call test_norm2_s_11d(error)
        if (error) return
        call test_norm2_s_12d(error)
        if (error) return
        call test_norm2_s_13d(error)
        if (error) return
        call test_norm2_s_14d(error)
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
        call test_norm_d_7d(error)
        if (error) return
        call test_norm_d_8d(error)
        if (error) return
        call test_norm_d_9d(error)
        if (error) return
        call test_norm_d_10d(error)
        if (error) return
        call test_norm_d_11d(error)
        if (error) return
        call test_norm_d_12d(error)
        if (error) return
        call test_norm_d_13d(error)
        if (error) return
        call test_norm_d_14d(error)
        if (error) return
        call test_norm2_d_2d(error)
        if (error) return
        call test_norm2_d_3d(error)
        if (error) return
        call test_norm2_d_4d(error)
        if (error) return
        call test_norm2_d_5d(error)
        if (error) return
        call test_norm2_d_6d(error)
        if (error) return
        call test_norm2_d_7d(error)
        if (error) return
        call test_norm2_d_8d(error)
        if (error) return
        call test_norm2_d_9d(error)
        if (error) return
        call test_norm2_d_10d(error)
        if (error) return
        call test_norm2_d_11d(error)
        if (error) return
        call test_norm2_d_12d(error)
        if (error) return
        call test_norm2_d_13d(error)
        if (error) return
        call test_norm2_d_14d(error)
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
        call test_norm_q_7d(error)
        if (error) return
        call test_norm_q_8d(error)
        if (error) return
        call test_norm_q_9d(error)
        if (error) return
        call test_norm_q_10d(error)
        if (error) return
        call test_norm_q_11d(error)
        if (error) return
        call test_norm_q_12d(error)
        if (error) return
        call test_norm_q_13d(error)
        if (error) return
        call test_norm_q_14d(error)
        if (error) return
        call test_norm2_q_2d(error)
        if (error) return
        call test_norm2_q_3d(error)
        if (error) return
        call test_norm2_q_4d(error)
        if (error) return
        call test_norm2_q_5d(error)
        if (error) return
        call test_norm2_q_6d(error)
        if (error) return
        call test_norm2_q_7d(error)
        if (error) return
        call test_norm2_q_8d(error)
        if (error) return
        call test_norm2_q_9d(error)
        if (error) return
        call test_norm2_q_10d(error)
        if (error) return
        call test_norm2_q_11d(error)
        if (error) return
        call test_norm2_q_12d(error)
        if (error) return
        call test_norm2_q_13d(error)
        if (error) return
        call test_norm2_q_14d(error)
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
        call test_norm_c_7d(error)
        if (error) return
        call test_norm_c_8d(error)
        if (error) return
        call test_norm_c_9d(error)
        if (error) return
        call test_norm_c_10d(error)
        if (error) return
        call test_norm_c_11d(error)
        if (error) return
        call test_norm_c_12d(error)
        if (error) return
        call test_norm_c_13d(error)
        if (error) return
        call test_norm_c_14d(error)
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
        call test_norm_z_7d(error)
        if (error) return
        call test_norm_z_8d(error)
        if (error) return
        call test_norm_z_9d(error)
        if (error) return
        call test_norm_z_10d(error)
        if (error) return
        call test_norm_z_11d(error)
        if (error) return
        call test_norm_z_12d(error)
        if (error) return
        call test_norm_z_13d(error)
        if (error) return
        call test_norm_z_14d(error)
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
        call test_norm_w_7d(error)
        if (error) return
        call test_norm_w_8d(error)
        if (error) return
        call test_norm_w_9d(error)
        if (error) return
        call test_norm_w_10d(error)
        if (error) return
        call test_norm_w_11d(error)
        if (error) return
        call test_norm_w_12d(error)
        if (error) return
        call test_norm_w_13d(error)
        if (error) return
        call test_norm_w_14d(error)
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
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**7
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_7d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**8
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_8d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**9
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_9d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**10
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_10d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**11
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_11d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**12
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_12d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**13
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_13d
    
    !> Test several norms with different dimensions
    subroutine test_norm_s_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**14
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_s_14d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_2d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_3d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_4d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_5d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_6d

    subroutine test_norm2_s_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 7
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_7d

    subroutine test_norm2_s_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 8
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_8d

    subroutine test_norm2_s_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 9
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_9d

    subroutine test_norm2_s_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 10
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_10d

    subroutine test_norm2_s_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 11
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_11d

    subroutine test_norm2_s_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 12
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_12d

    subroutine test_norm2_s_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 13
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_13d

    subroutine test_norm2_s_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 14
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        real(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_sp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_sp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_s_14d

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
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**7
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_7d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**8
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_8d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**9
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_9d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**10
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_10d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**11
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_11d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**12
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_12d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**13
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_13d
    
    !> Test several norms with different dimensions
    subroutine test_norm_d_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**14
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_d_14d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_2d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_3d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_4d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_5d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_6d

    subroutine test_norm2_d_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 7
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_7d

    subroutine test_norm2_d_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 8
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_8d

    subroutine test_norm2_d_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 9
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_9d

    subroutine test_norm2_d_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 10
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_10d

    subroutine test_norm2_d_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 11
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_11d

    subroutine test_norm2_d_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 12
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_12d

    subroutine test_norm2_d_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 13
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_13d

    subroutine test_norm2_d_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 14
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        real(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_dp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_dp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_d_14d

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
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**7
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_7d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**8
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_8d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**9
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_9d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**10
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_10d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**11
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_11d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**12
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_12d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**13
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_13d
    
    !> Test several norms with different dimensions
    subroutine test_norm_q_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**14
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_q_14d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_2d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_3d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_4d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_5d

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
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_6d

    subroutine test_norm2_q_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 7
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_7d

    subroutine test_norm2_q_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 8
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_8d

    subroutine test_norm2_q_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 9
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_9d

    subroutine test_norm2_q_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 10
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_10d

    subroutine test_norm2_q_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 11
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_11d

    subroutine test_norm2_q_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 12
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_12d

    subroutine test_norm2_q_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 13
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_13d

    subroutine test_norm2_q_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,dim
        integer(ilp),parameter :: ndim = 14
        integer(ilp),parameter :: n = 2_ilp**ndim
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        real(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        intrinsic :: norm2
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
        ! Init as a range,but with small elements such that all power norms will
        ! never overflow, even in single precision
        a = [(0.01_qp*(j - n/2_ilp),j=1_ilp,n)]
        b = reshape(a,shape(b))
        
        ! Test some norms
        error = .not. abs(norm(a,2) - norm2(a)) < tol*norm(a,2)
        print *, 'a ',norm(a,2),norm2(a)
        if (error) return
        
        ! Infinity norms
        error = .not. abs(norm(b,2) - norm2(b)) < tol*norm(b,2)
        print *, 'B ',norm(b,2),norm2(b)
        if (error) return
        
        ! Test norm as collapsed around dimension
        do dim = 1,ndim
            
            error = .not. all(abs(norm(b,2,dim) - norm2(b,dim)) < tol*max(1.0_qp,norm(b,2,dim)))
            if (error) then
            print *, 'dim=',dim
            print *, ' normA=',norm(b,2,dim)
            print *, ' normB=',norm2(b,dim)
            if (error) return
            end if
            
        end do
        
    end subroutine test_norm2_q_14d

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
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**7
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_7d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**8
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_8d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**9
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_9d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**10
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_10d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**11
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_11d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**12
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_12d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**13
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_13d
    
    !> Test several norms with different dimensions
    subroutine test_norm_c_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**14
        real(sp),parameter :: tol = 10*sqrt(epsilon(0.0_sp))
        complex(sp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_c_14d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals

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
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**7
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_7d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**8
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_8d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**9
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_9d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**10
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_10d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**11
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_11d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**12
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_12d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**13
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_13d
    
    !> Test several norms with different dimensions
    subroutine test_norm_z_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**14
        real(dp),parameter :: tol = 10*sqrt(epsilon(0.0_dp))
        complex(dp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_z_14d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals

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
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_7d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**7
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_7d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_8d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**8
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_8d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_9d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**9
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_9d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_10d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**10
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_10d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_11d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**11
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_11d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_12d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**12
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_12d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_13d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**13
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_13d
    
    !> Test several norms with different dimensions
    subroutine test_norm_w_14d(error)
        logical,intent(out) :: error

        integer(ilp) :: j,order
        integer(ilp),parameter :: n = 2_ilp**14
        real(qp),parameter :: tol = 10*sqrt(epsilon(0.0_qp))
        complex(qp),allocatable :: a(:),b(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        
        allocate (a(n),b(2,2,2,2,2,2,2,2,2,2,2,2,2,2))
        
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
        
    end subroutine test_norm_w_14d

    !> Test Euclidean norm; compare with Fortran intrinsic norm2 for reals

end module test_linalg_norms

