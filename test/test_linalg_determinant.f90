! Test matrix determinant
module test_linalg_determinant
    use linear_algebra

    implicit none(type,external)

    contains

    !> Matrix inversion tests
    subroutine test_matrix_determinant(error)
        logical,intent(out) :: error

        real :: t0,t1

        call cpu_time(t0)

        call test_s_eye_determinant(error)
        if (error) return

        call test_s_eye_multiple(error)
        if (error) return
        call test_d_eye_determinant(error)
        if (error) return

        call test_d_eye_multiple(error)
        if (error) return
        call test_q_eye_determinant(error)
        if (error) return

        call test_q_eye_multiple(error)
        if (error) return
        call test_c_eye_determinant(error)
        if (error) return

        call test_c_eye_multiple(error)
        if (error) return
        call test_z_eye_determinant(error)
        if (error) return

        call test_z_eye_multiple(error)
        if (error) return
        call test_w_eye_determinant(error)
        if (error) return

        call test_w_eye_multiple(error)
        if (error) return

        call test_c_complex_determinant(error)
        if (error) return
        call test_z_complex_determinant(error)
        if (error) return
        call test_w_complex_determinant(error)
        if (error) return

        call cpu_time(t1)

        print 1,1000*(t1 - t0),merge('SUCCESS','ERROR  ',.not. error)

1       format('Determinant tests completed in ',f9.4,' milliseconds, result=',a)

    end subroutine test_matrix_determinant

    !> Determinant of identity matrix
    subroutine test_s_eye_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i
        integer(ilp),parameter :: n = 128_ilp

        real(sp) :: a(n,n),deta

        a = eye(n)

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_sp) < tiny(0.0_sp)
        if (error) return

    end subroutine test_s_eye_determinant

    !> Determinant of identity matrix multiplier
    subroutine test_s_eye_multiple(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(sp),parameter :: coef = 0.0001_sp
        integer(ilp) :: i,j
        real(sp) :: a(n,n),deta

        !> Multiply eye by a very small number
        a = eye(n)
        do concurrent(i=1:n)
          a(i,i) = coef
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < max(tiny(0.0_sp),epsilon(0.0_sp)*coef**n)
        if (error) return

    end subroutine test_s_eye_multiple

    subroutine test_d_eye_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i
        integer(ilp),parameter :: n = 128_ilp

        real(dp) :: a(n,n),deta

        a = eye(n)

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_dp) < tiny(0.0_dp)
        if (error) return

    end subroutine test_d_eye_determinant

    !> Determinant of identity matrix multiplier
    subroutine test_d_eye_multiple(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(dp),parameter :: coef = 0.0001_dp
        integer(ilp) :: i,j
        real(dp) :: a(n,n),deta

        !> Multiply eye by a very small number
        a = eye(n)
        do concurrent(i=1:n)
          a(i,i) = coef
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < max(tiny(0.0_dp),epsilon(0.0_dp)*coef**n)
        if (error) return

    end subroutine test_d_eye_multiple

    subroutine test_q_eye_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i
        integer(ilp),parameter :: n = 128_ilp

        real(qp) :: a(n,n),deta

        a = eye(n)

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_qp) < tiny(0.0_qp)
        if (error) return

    end subroutine test_q_eye_determinant

    !> Determinant of identity matrix multiplier
    subroutine test_q_eye_multiple(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(qp),parameter :: coef = 0.0001_qp
        integer(ilp) :: i,j
        real(qp) :: a(n,n),deta

        !> Multiply eye by a very small number
        a = eye(n)
        do concurrent(i=1:n)
          a(i,i) = coef
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < max(tiny(0.0_qp),epsilon(0.0_qp)*coef**n)
        if (error) return

    end subroutine test_q_eye_multiple

    subroutine test_c_eye_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i
        integer(ilp),parameter :: n = 128_ilp

        complex(sp) :: a(n,n),deta

        a = eye(n)

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_sp) < tiny(0.0_sp)
        if (error) return

    end subroutine test_c_eye_determinant

    !> Determinant of identity matrix multiplier
    subroutine test_c_eye_multiple(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(sp),parameter :: coef = 0.0001_sp
        integer(ilp) :: i,j
        complex(sp) :: a(n,n),deta

        !> Multiply eye by a very small number
        a = eye(n)
        do concurrent(i=1:n)
          a(i,i) = coef
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < max(tiny(0.0_sp),epsilon(0.0_sp)*coef**n)
        if (error) return

    end subroutine test_c_eye_multiple

    subroutine test_z_eye_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i
        integer(ilp),parameter :: n = 128_ilp

        complex(dp) :: a(n,n),deta

        a = eye(n)

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_dp) < tiny(0.0_dp)
        if (error) return

    end subroutine test_z_eye_determinant

    !> Determinant of identity matrix multiplier
    subroutine test_z_eye_multiple(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(dp),parameter :: coef = 0.0001_dp
        integer(ilp) :: i,j
        complex(dp) :: a(n,n),deta

        !> Multiply eye by a very small number
        a = eye(n)
        do concurrent(i=1:n)
          a(i,i) = coef
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < max(tiny(0.0_dp),epsilon(0.0_dp)*coef**n)
        if (error) return

    end subroutine test_z_eye_multiple

    subroutine test_w_eye_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i
        integer(ilp),parameter :: n = 128_ilp

        complex(qp) :: a(n,n),deta

        a = eye(n)

        !> Determinant function
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - 1.0_qp) < tiny(0.0_qp)
        if (error) return

    end subroutine test_w_eye_determinant

    !> Determinant of identity matrix multiplier
    subroutine test_w_eye_multiple(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 10_ilp
        real(qp),parameter :: coef = 0.0001_qp
        integer(ilp) :: i,j
        complex(qp) :: a(n,n),deta

        !> Multiply eye by a very small number
        a = eye(n)
        do concurrent(i=1:n)
          a(i,i) = coef
        end do

        !> Determinant: small, but a is not singular, because it is a multiple of the identity.
        deta = det(a,err=state)
        error = state%error() .or. .not. abs(deta - coef**n) < max(tiny(0.0_qp),epsilon(0.0_qp)*coef**n)
        if (error) return

    end subroutine test_w_eye_multiple

    !> Determinant of complex identity matrix
    subroutine test_c_complex_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,n
        integer(ilp),parameter :: nmax = 10_ilp

        complex(sp),parameter :: res(nmax) = [complex(sp) ::(1,1), (0,2), (-2,2), (-4,0), (-4,-4), &
                                                  (0,-8), (8,-8), (16,0), (16,16), (0,32)]

        complex(sp),allocatable :: a(:,:)
        complex(sp) :: deta(nmax)

        !> Test determinant for all sizes, 1:nmax
        matrix_size: do n = 1,nmax

           ! Put 1+i on each diagonal element
           a = eye(n)
           do concurrent(i=1:n)
             a(i,i) = (1.0_sp,1.0_sp)
           end do

           ! Expected result
           deta(n) = det(a,err=state)

           deallocate (a)
           if (state%error()) exit matrix_size

        end do matrix_size

        error = state%error() .or. any(.not. abs(res - deta) <= tiny(0.0_sp))

    end subroutine test_c_complex_determinant

    subroutine test_z_complex_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,n
        integer(ilp),parameter :: nmax = 10_ilp

        complex(dp),parameter :: res(nmax) = [complex(dp) ::(1,1), (0,2), (-2,2), (-4,0), (-4,-4), &
                                                  (0,-8), (8,-8), (16,0), (16,16), (0,32)]

        complex(dp),allocatable :: a(:,:)
        complex(dp) :: deta(nmax)

        !> Test determinant for all sizes, 1:nmax
        matrix_size: do n = 1,nmax

           ! Put 1+i on each diagonal element
           a = eye(n)
           do concurrent(i=1:n)
             a(i,i) = (1.0_dp,1.0_dp)
           end do

           ! Expected result
           deta(n) = det(a,err=state)

           deallocate (a)
           if (state%error()) exit matrix_size

        end do matrix_size

        error = state%error() .or. any(.not. abs(res - deta) <= tiny(0.0_dp))

    end subroutine test_z_complex_determinant

    subroutine test_w_complex_determinant(error)
        logical,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,n
        integer(ilp),parameter :: nmax = 10_ilp

        complex(qp),parameter :: res(nmax) = [complex(qp) ::(1,1), (0,2), (-2,2), (-4,0), (-4,-4), &
                                                  (0,-8), (8,-8), (16,0), (16,16), (0,32)]

        complex(qp),allocatable :: a(:,:)
        complex(qp) :: deta(nmax)

        !> Test determinant for all sizes, 1:nmax
        matrix_size: do n = 1,nmax

           ! Put 1+i on each diagonal element
           a = eye(n)
           do concurrent(i=1:n)
             a(i,i) = (1.0_qp,1.0_qp)
           end do

           ! Expected result
           deta(n) = det(a,err=state)

           deallocate (a)
           if (state%error()) exit matrix_size

        end do matrix_size

        error = state%error() .or. any(.not. abs(res - deta) <= tiny(0.0_qp))

    end subroutine test_w_complex_determinant

end module test_linalg_determinant

