! Test inverse matrix operator
module test_linalg_inverse
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants
    use linear_algebra,only:inv,invert,operator(.inv.),eye
    use la_state_type,only:la_state,LINALG_ERROR

    implicit none(type,external)
    private
    
    public :: test_inverse_matrix

    contains

    !> Matrix inversion tests
    subroutine test_inverse_matrix(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)
        
        allocate (tests(0))

        call add_test(tests,new_unittest("s_eye_inverse",test_s_eye_inverse))
        call add_test(tests,new_unittest("s_singular_inverse",test_s_singular_inverse))
        call add_test(tests,new_unittest("s_random_spd_inverse",test_s_random_spd_inverse))
        call add_test(tests,new_unittest("d_eye_inverse",test_d_eye_inverse))
        call add_test(tests,new_unittest("d_singular_inverse",test_d_singular_inverse))
        call add_test(tests,new_unittest("d_random_spd_inverse",test_d_random_spd_inverse))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("x_eye_inverse",test_x_eye_inverse))
        call add_test(tests,new_unittest("x_singular_inverse",test_x_singular_inverse))
        call add_test(tests,new_unittest("x_random_spd_inverse",test_x_random_spd_inverse))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("q_eye_inverse",test_q_eye_inverse))
        call add_test(tests,new_unittest("q_singular_inverse",test_q_singular_inverse))
        call add_test(tests,new_unittest("q_random_spd_inverse",test_q_random_spd_inverse))
#endif
        call add_test(tests,new_unittest("c_eye_inverse",test_c_eye_inverse))
        call add_test(tests,new_unittest("c_singular_inverse",test_c_singular_inverse))
        call add_test(tests,new_unittest("c_random_spd_inverse",test_c_random_spd_inverse))
        call add_test(tests,new_unittest("z_eye_inverse",test_z_eye_inverse))
        call add_test(tests,new_unittest("z_singular_inverse",test_z_singular_inverse))
        call add_test(tests,new_unittest("z_random_spd_inverse",test_z_random_spd_inverse))
#ifdef LA_WITH_XDP
        call add_test(tests,new_unittest("y_eye_inverse",test_y_eye_inverse))
        call add_test(tests,new_unittest("y_singular_inverse",test_y_singular_inverse))
        call add_test(tests,new_unittest("y_random_spd_inverse",test_y_random_spd_inverse))
#endif
#ifdef LA_WITH_QP
        call add_test(tests,new_unittest("w_eye_inverse",test_w_eye_inverse))
        call add_test(tests,new_unittest("w_singular_inverse",test_w_singular_inverse))
        call add_test(tests,new_unittest("w_random_spd_inverse",test_w_random_spd_inverse))
#endif

    end subroutine test_inverse_matrix

    !> Invert real identity matrix
    subroutine test_s_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 25_ilp
        real(sp) :: a(n,n),inva(n,n)

        a = real(eye(n),sp)

        !> Inverse function
        inva = inv(a,err=state)
        call check(error,state%ok(),'inverse_s_eye (function): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_sp)),'inverse_s_eye (function): data converged')
        if (allocated(error)) return
        
        !> Inverse subroutine in-place
        call invert(a,err=state)

        call check(error,state%ok(),'inverse_s_eye (in-place): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_sp)),'inverse_s_eye (in-place): data converged')
        if (allocated(error)) return

    end subroutine test_s_eye_inverse

    !> Invert singular matrix
    subroutine test_s_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        real(sp) :: a(n,n)

        a = real(eye(n),sp)
        
        !> Make rank-deficient
        a(12,12) = 0
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular real(sp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_s_singular_inverse
    
    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_s(n) result(A)
        integer(ilp),intent(in) :: n
        real(sp) :: A(n,n)
        
        real(sp),parameter :: half = 0.5_sp
        
        !> Initialize with randoms
        call random_number(A)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        A = real(A + n*eye(n),sp)
        
    end function random_spd_matrix_s

    !> Test random symmetric positive-definite matrix
    subroutine test_s_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Solution tolerance
        real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))

        !> Local variables
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        real(sp) :: A(n,n),Am1(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_s(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random SPD matrix (sp): '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(matmul(Am1,A) - eye(n)) < tol),'random SPD matrix (sp): accuracy test')
        if (allocated(error)) return
    
    end subroutine test_s_random_spd_inverse
    
    !> Invert real identity matrix
    subroutine test_d_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 25_ilp
        real(dp) :: a(n,n),inva(n,n)

        a = real(eye(n),dp)

        !> Inverse function
        inva = inv(a,err=state)
        call check(error,state%ok(),'inverse_d_eye (function): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_dp)),'inverse_d_eye (function): data converged')
        if (allocated(error)) return
        
        !> Inverse subroutine in-place
        call invert(a,err=state)

        call check(error,state%ok(),'inverse_d_eye (in-place): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_dp)),'inverse_d_eye (in-place): data converged')
        if (allocated(error)) return

    end subroutine test_d_eye_inverse

    !> Invert singular matrix
    subroutine test_d_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        real(dp) :: a(n,n)

        a = real(eye(n),dp)
        
        !> Make rank-deficient
        a(12,12) = 0
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular real(dp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_d_singular_inverse
    
    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_d(n) result(A)
        integer(ilp),intent(in) :: n
        real(dp) :: A(n,n)
        
        real(dp),parameter :: half = 0.5_dp
        
        !> Initialize with randoms
        call random_number(A)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        A = real(A + n*eye(n),dp)
        
    end function random_spd_matrix_d

    !> Test random symmetric positive-definite matrix
    subroutine test_d_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Solution tolerance
        real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))

        !> Local variables
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        real(dp) :: A(n,n),Am1(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_d(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random SPD matrix (dp): '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(matmul(Am1,A) - eye(n)) < tol),'random SPD matrix (dp): accuracy test')
        if (allocated(error)) return
    
    end subroutine test_d_random_spd_inverse
    
#ifdef LA_WITH_XDP
    !> Invert real identity matrix
    subroutine test_x_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 25_ilp
        real(xdp) :: a(n,n),inva(n,n)

        a = real(eye(n),xdp)

        !> Inverse function
        inva = inv(a,err=state)
        call check(error,state%ok(),'inverse_x_eye (function): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_xdp)),'inverse_x_eye (function): data converged')
        if (allocated(error)) return
        
        !> Inverse subroutine in-place
        call invert(a,err=state)

        call check(error,state%ok(),'inverse_x_eye (in-place): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_xdp)),'inverse_x_eye (in-place): data converged')
        if (allocated(error)) return

    end subroutine test_x_eye_inverse

    !> Invert singular matrix
    subroutine test_x_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        real(xdp) :: a(n,n)

        a = real(eye(n),xdp)
        
        !> Make rank-deficient
        a(12,12) = 0
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular real(xdp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_x_singular_inverse
    
    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_x(n) result(A)
        integer(ilp),intent(in) :: n
        real(xdp) :: A(n,n)
        
        real(xdp),parameter :: half = 0.5_xdp
        
        !> Initialize with randoms
        call random_number(A)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        A = real(A + n*eye(n),xdp)
        
    end function random_spd_matrix_x

    !> Test random symmetric positive-definite matrix
    subroutine test_x_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Solution tolerance
        real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))

        !> Local variables
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        real(xdp) :: A(n,n),Am1(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_x(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random SPD matrix (xdp): '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(matmul(Am1,A) - eye(n)) < tol),'random SPD matrix (xdp): accuracy test')
        if (allocated(error)) return
    
    end subroutine test_x_random_spd_inverse
#endif
    
#ifdef LA_WITH_QP
    !> Invert real identity matrix
    subroutine test_q_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp),parameter :: n = 25_ilp
        real(qp) :: a(n,n),inva(n,n)

        a = real(eye(n),qp)

        !> Inverse function
        inva = inv(a,err=state)
        call check(error,state%ok(),'inverse_q_eye (function): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_qp)),'inverse_q_eye (function): data converged')
        if (allocated(error)) return
        
        !> Inverse subroutine in-place
        call invert(a,err=state)

        call check(error,state%ok(),'inverse_q_eye (in-place): '//state%print())
        if (allocated(error)) return
        
        call check(error,all(abs(a - inva) < epsilon(0.0_qp)),'inverse_q_eye (in-place): data converged')
        if (allocated(error)) return

    end subroutine test_q_eye_inverse

    !> Invert singular matrix
    subroutine test_q_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        real(qp) :: a(n,n)

        a = real(eye(n),qp)
        
        !> Make rank-deficient
        a(12,12) = 0
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular real(qp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_q_singular_inverse
    
    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_q(n) result(A)
        integer(ilp),intent(in) :: n
        real(qp) :: A(n,n)
        
        real(qp),parameter :: half = 0.5_qp
        
        !> Initialize with randoms
        call random_number(A)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        A = real(A + n*eye(n),qp)
        
    end function random_spd_matrix_q

    !> Test random symmetric positive-definite matrix
    subroutine test_q_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Solution tolerance
        real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))

        !> Local variables
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        real(qp) :: A(n,n),Am1(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_q(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random SPD matrix (qp): '//state%print())
        if (allocated(error)) return

        call check(error,all(abs(matmul(Am1,A) - eye(n)) < tol),'random SPD matrix (qp): accuracy test')
        if (allocated(error)) return
    
    end subroutine test_q_random_spd_inverse
#endif
    
    !> Invert complex identity matrix
    subroutine test_c_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 25_ilp

        complex(sp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_sp,1.0_sp), (0.0_sp,0.0_sp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        call check(error,state%ok(),'inverse_c_eye (function): '//state%print())
        if (allocated(error)) return
        
        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_c_eye (function): data converged')
        if (allocated(error)) return

        !> Inverse subroutine
        call invert(copya,err=state)

        call check(error,state%ok(),'inverse_c_eye (subroutine): '//state%print())
        if (allocated(error)) return

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_c_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(sp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < epsilon(0.0_sp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_sp,0.0_sp)) < epsilon(0.0_sp)
               end if
           end function is_diagonal_inverse

    end subroutine test_c_eye_inverse

    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_c(n) result(A)
        integer(ilp),intent(in) :: n
        complex(sp) :: A(n,n)
        
        complex(sp),parameter :: half = (0.5_sp,0.0_sp)
        real(sp) :: reA(n,n),imA(n,n)
        integer(ilp) :: i
        
        !> Initialize with randoms
        call random_number(reA)
        call random_number(imA)
        
        A = cmplx(reA,imA,kind=sp)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        forall (i=1:n) A(i,i) = A(i,i) + n*(1.0_sp,0.0_sp)
        
    end function random_spd_matrix_c

    !> Test random symmetric positive-definite matrix
    subroutine test_c_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Local variables
        integer(ilp) :: failed,i,j
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        complex(sp) :: A(n,n),Am1(n,n),AA(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_c(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random complex SPD matrix (sp): '//state%print())
        if (allocated(error)) return

        failed = 0
        AA = matmul(A,Am1)
        do i = 1,n
            do j = 1,n
                if (.not. is_complex_inverse(AA(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_c_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_complex_inverse(aij,i,j)
               complex(sp),intent(in) :: aij
               integer(ilp),intent(in) :: i,j
               real(sp),parameter :: tol = sqrt(epsilon(0.0_sp))
               if (i /= j) then
                  is_complex_inverse = abs(aij) < tol
               else
                  ! Product should return the real identity
                  is_complex_inverse = abs(aij - (1.0_sp,0.0_sp)) < tol
               end if
           end function is_complex_inverse
    
    end subroutine test_c_random_spd_inverse

    !> Invert singular matrix
    subroutine test_c_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        complex(sp) :: a(n,n)

        a = (0.0_sp,0.0_sp)
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular complex(sp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_c_singular_inverse

    subroutine test_z_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 25_ilp

        complex(dp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_dp,1.0_dp), (0.0_dp,0.0_dp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        call check(error,state%ok(),'inverse_z_eye (function): '//state%print())
        if (allocated(error)) return
        
        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_z_eye (function): data converged')
        if (allocated(error)) return

        !> Inverse subroutine
        call invert(copya,err=state)

        call check(error,state%ok(),'inverse_z_eye (subroutine): '//state%print())
        if (allocated(error)) return

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_z_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(dp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < epsilon(0.0_dp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_dp,0.0_dp)) < epsilon(0.0_dp)
               end if
           end function is_diagonal_inverse

    end subroutine test_z_eye_inverse

    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_z(n) result(A)
        integer(ilp),intent(in) :: n
        complex(dp) :: A(n,n)
        
        complex(dp),parameter :: half = (0.5_dp,0.0_dp)
        real(dp) :: reA(n,n),imA(n,n)
        integer(ilp) :: i
        
        !> Initialize with randoms
        call random_number(reA)
        call random_number(imA)
        
        A = cmplx(reA,imA,kind=dp)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        forall (i=1:n) A(i,i) = A(i,i) + n*(1.0_dp,0.0_dp)
        
    end function random_spd_matrix_z

    !> Test random symmetric positive-definite matrix
    subroutine test_z_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Local variables
        integer(ilp) :: failed,i,j
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        complex(dp) :: A(n,n),Am1(n,n),AA(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_z(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random complex SPD matrix (dp): '//state%print())
        if (allocated(error)) return

        failed = 0
        AA = matmul(A,Am1)
        do i = 1,n
            do j = 1,n
                if (.not. is_complex_inverse(AA(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_z_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_complex_inverse(aij,i,j)
               complex(dp),intent(in) :: aij
               integer(ilp),intent(in) :: i,j
               real(dp),parameter :: tol = sqrt(epsilon(0.0_dp))
               if (i /= j) then
                  is_complex_inverse = abs(aij) < tol
               else
                  ! Product should return the real identity
                  is_complex_inverse = abs(aij - (1.0_dp,0.0_dp)) < tol
               end if
           end function is_complex_inverse
    
    end subroutine test_z_random_spd_inverse

    !> Invert singular matrix
    subroutine test_z_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        complex(dp) :: a(n,n)

        a = (0.0_dp,0.0_dp)
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular complex(dp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_z_singular_inverse

#ifdef LA_WITH_XDP
    subroutine test_y_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 25_ilp

        complex(xdp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_xdp,1.0_xdp), (0.0_xdp,0.0_xdp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        call check(error,state%ok(),'inverse_y_eye (function): '//state%print())
        if (allocated(error)) return
        
        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_y_eye (function): data converged')
        if (allocated(error)) return

        !> Inverse subroutine
        call invert(copya,err=state)

        call check(error,state%ok(),'inverse_y_eye (subroutine): '//state%print())
        if (allocated(error)) return

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_y_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(xdp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < epsilon(0.0_xdp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_xdp,0.0_xdp)) < epsilon(0.0_xdp)
               end if
           end function is_diagonal_inverse

    end subroutine test_y_eye_inverse

    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_y(n) result(A)
        integer(ilp),intent(in) :: n
        complex(xdp) :: A(n,n)
        
        complex(xdp),parameter :: half = (0.5_xdp,0.0_xdp)
        real(xdp) :: reA(n,n),imA(n,n)
        integer(ilp) :: i
        
        !> Initialize with randoms
        call random_number(reA)
        call random_number(imA)
        
        A = cmplx(reA,imA,kind=xdp)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        forall (i=1:n) A(i,i) = A(i,i) + n*(1.0_xdp,0.0_xdp)
        
    end function random_spd_matrix_y

    !> Test random symmetric positive-definite matrix
    subroutine test_y_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Local variables
        integer(ilp) :: failed,i,j
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        complex(xdp) :: A(n,n),Am1(n,n),AA(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_y(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random complex SPD matrix (xdp): '//state%print())
        if (allocated(error)) return

        failed = 0
        AA = matmul(A,Am1)
        do i = 1,n
            do j = 1,n
                if (.not. is_complex_inverse(AA(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_y_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_complex_inverse(aij,i,j)
               complex(xdp),intent(in) :: aij
               integer(ilp),intent(in) :: i,j
               real(xdp),parameter :: tol = sqrt(epsilon(0.0_xdp))
               if (i /= j) then
                  is_complex_inverse = abs(aij) < tol
               else
                  ! Product should return the real identity
                  is_complex_inverse = abs(aij - (1.0_xdp,0.0_xdp)) < tol
               end if
           end function is_complex_inverse
    
    end subroutine test_y_random_spd_inverse

    !> Invert singular matrix
    subroutine test_y_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        complex(xdp) :: a(n,n)

        a = (0.0_xdp,0.0_xdp)
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular complex(xdp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_y_singular_inverse
#endif

#ifdef LA_WITH_QP
    subroutine test_w_eye_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: state

        integer(ilp) :: i,j,failed
        integer(ilp),parameter :: n = 25_ilp

        complex(qp) :: a(n,n),copya(n,n),inva(n,n)

        do concurrent(i=1:n,j=1:n)
          a(i,j) = merge((1.0_qp,1.0_qp), (0.0_qp,0.0_qp),i == j)
        end do
        copya = a

        !> The inverse of a complex diagonal matrix has conjg(z_ii)/abs(z_ii)^2 on the diagonal
        inva = inv(a,err=state)

        call check(error,state%ok(),'inverse_w_eye (function): '//state%print())
        if (allocated(error)) return
        
        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),inva(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_w_eye (function): data converged')
        if (allocated(error)) return

        !> Inverse subroutine
        call invert(copya,err=state)

        call check(error,state%ok(),'inverse_w_eye (subroutine): '//state%print())
        if (allocated(error)) return

        failed = 0
        do i = 1,n
            do j = 1,n
                if (.not. is_diagonal_inverse(a(i,j),copya(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_w_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_diagonal_inverse(aij,invaij,i,j)
               complex(qp),intent(in) :: aij,invaij
               integer(ilp),intent(in) :: i,j
               if (i /= j) then
                  is_diagonal_inverse = max(abs(aij),abs(invaij)) < epsilon(0.0_qp)
               else
                  ! Product should return the real identity
                  is_diagonal_inverse = abs(aij*invaij - (1.0_qp,0.0_qp)) < epsilon(0.0_qp)
               end if
           end function is_diagonal_inverse

    end subroutine test_w_eye_inverse

    !> Create a random symmetric positive definite matrix
    function random_spd_matrix_w(n) result(A)
        integer(ilp),intent(in) :: n
        complex(qp) :: A(n,n)
        
        complex(qp),parameter :: half = (0.5_qp,0.0_qp)
        real(qp) :: reA(n,n),imA(n,n)
        integer(ilp) :: i
        
        !> Initialize with randoms
        call random_number(reA)
        call random_number(imA)
        
        A = cmplx(reA,imA,kind=qp)
        
        !> Make symmetric
        A = half*(A + transpose(A))
        
        !> Add diagonally dominant part
        forall (i=1:n) A(i,i) = A(i,i) + n*(1.0_qp,0.0_qp)
        
    end function random_spd_matrix_w

    !> Test random symmetric positive-definite matrix
    subroutine test_w_random_spd_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        !> Local variables
        integer(ilp) :: failed,i,j
        integer(ilp),parameter :: n = 5_ilp
        type(la_state) :: state
        complex(qp) :: A(n,n),Am1(n,n),AA(n,n)
        
        !> Generate random SPD matrix
        A = random_spd_matrix_w(n)

        !> Invert matrix
        Am1 = inv(A,err=state)
        
        !> Check result
        call check(error,state%ok(),'random complex SPD matrix (qp): '//state%print())
        if (allocated(error)) return

        failed = 0
        AA = matmul(A,Am1)
        do i = 1,n
            do j = 1,n
                if (.not. is_complex_inverse(AA(i,j),i,j)) failed = failed + 1
            end do
        end do

        call check(error,failed == 0,'inverse_w_eye (subroutine): data converged')
        if (allocated(error)) return

        contains

           elemental logical function is_complex_inverse(aij,i,j)
               complex(qp),intent(in) :: aij
               integer(ilp),intent(in) :: i,j
               real(qp),parameter :: tol = sqrt(epsilon(0.0_qp))
               if (i /= j) then
                  is_complex_inverse = abs(aij) < tol
               else
                  ! Product should return the real identity
                  is_complex_inverse = abs(aij - (1.0_qp,0.0_qp)) < tol
               end if
           end function is_complex_inverse
    
    end subroutine test_w_random_spd_inverse

    !> Invert singular matrix
    subroutine test_w_singular_inverse(error)
        type(error_type),allocatable,intent(out) :: error

        type(la_state) :: err

        integer(ilp),parameter :: n = 25_ilp
        complex(qp) :: a(n,n)

        a = (0.0_qp,0.0_qp)
        
        !> Inverse function
        call invert(a,err=err)
        call check(error,err%state == LINALG_ERROR,'singular complex(qp) inverse returned '//err%print())
        if (allocated(error)) return
        
    end subroutine test_w_singular_inverse
#endif

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type),allocatable,intent(inout) :: tests(:)
        type(unittest_type),intent(in) :: new_test
        
        integer :: n
        type(unittest_type),allocatable :: new_tests(:)
        
        if (allocated(tests)) then
            n = size(tests)
        else
            n = 0
        end if
        
        allocate (new_tests(n + 1))
        if (n > 0) new_tests(1:n) = tests(1:n)
                 new_tests(1 + n) = new_test
        call move_alloc(from=new_tests,to=tests)
        
    end subroutine add_test

end module test_linalg_inverse

program test_inv
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_inverse,only:test_inverse_matrix
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_inverse",test_inverse_matrix) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_inv
