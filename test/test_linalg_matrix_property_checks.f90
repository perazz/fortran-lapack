
module test_linalg_matrix_property_checks
    use testdrive,only:new_unittest,unittest_type,error_type,check
    use la_constants,only:sp,dp,qp
    use linear_algebra,only:is_square,is_diagonal,is_symmetric, &
        is_skew_symmetric,is_hermitian,is_triangular,is_hessenberg

    implicit none

    real(sp),parameter :: sptol = 1000*epsilon(1._sp)
    real(dp),parameter :: dptol = 1000*epsilon(1._dp)
    real(qp),parameter :: qptol = 1000*epsilon(1._qp)

contains

    !> Collect all exported unit tests
    subroutine collect_linalg_matrix_property_checks(testsuite)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("is_square_rsp",test_is_square_rsp), &
            new_unittest("is_square_rdp",test_is_square_rdp), &
            new_unittest("is_square_rqp",test_is_square_rqp), &
            new_unittest("is_square_csp",test_is_square_csp), &
            new_unittest("is_square_cdp",test_is_square_cdp), &
            new_unittest("is_square_cqp",test_is_square_cqp), &
            new_unittest("is_diagonal_rsp",test_is_diagonal_rsp), &
            new_unittest("is_diagonal_rdp",test_is_diagonal_rdp), &
            new_unittest("is_diagonal_rqp",test_is_diagonal_rqp), &
            new_unittest("is_diagonal_csp",test_is_diagonal_csp), &
            new_unittest("is_diagonal_cdp",test_is_diagonal_cdp), &
            new_unittest("is_diagonal_cqp",test_is_diagonal_cqp), &
            new_unittest("is_symmetric_rsp",test_is_symmetric_rsp), &
            new_unittest("is_symmetric_rdp",test_is_symmetric_rdp), &
            new_unittest("is_symmetric_rqp",test_is_symmetric_rqp), &
            new_unittest("is_symmetric_csp",test_is_symmetric_csp), &
            new_unittest("is_symmetric_cdp",test_is_symmetric_cdp), &
            new_unittest("is_symmetric_cqp",test_is_symmetric_cqp), &
            new_unittest("is_skew_symmetric_rsp",test_is_skew_symmetric_rsp), &
            new_unittest("is_skew_symmetric_rdp",test_is_skew_symmetric_rdp), &
            new_unittest("is_skew_symmetric_rqp",test_is_skew_symmetric_rqp), &
            new_unittest("is_skew_symmetric_csp",test_is_skew_symmetric_csp), &
            new_unittest("is_skew_symmetric_cdp",test_is_skew_symmetric_cdp), &
            new_unittest("is_skew_symmetric_cqp",test_is_skew_symmetric_cqp), &
            new_unittest("is_hermitian_rsp",test_is_hermitian_rsp), &
            new_unittest("is_hermitian_rdp",test_is_hermitian_rdp), &
            new_unittest("is_hermitian_rqp",test_is_hermitian_rqp), &
            new_unittest("is_hermitian_csp",test_is_hermitian_csp), &
            new_unittest("is_hermitian_cdp",test_is_hermitian_cdp), &
            new_unittest("is_hermitian_cqp",test_is_hermitian_cqp), &
            new_unittest("is_triangular_rsp",test_is_triangular_rsp), &
            new_unittest("is_triangular_rdp",test_is_triangular_rdp), &
            new_unittest("is_triangular_rqp",test_is_triangular_rqp), &
            new_unittest("is_triangular_csp",test_is_triangular_csp), &
            new_unittest("is_triangular_cdp",test_is_triangular_cdp), &
            new_unittest("is_triangular_cqp",test_is_triangular_cqp), &
            new_unittest("is_hessenberg_rsp",test_is_hessenberg_rsp), &
            new_unittest("is_hessenberg_rdp",test_is_hessenberg_rdp), &
            new_unittest("is_hessenberg_rqp",test_is_hessenberg_rqp), &
            new_unittest("is_hessenberg_csp",test_is_hessenberg_csp), &
            new_unittest("is_hessenberg_cdp",test_is_hessenberg_cdp), &
            new_unittest("is_hessenberg_cqp",test_is_hessenberg_cqp) &
            ]
    end subroutine collect_linalg_matrix_property_checks

    !is_square
    subroutine test_is_square_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true(2,2),A_false(2,3)
        A_true = reshape([1.,2.,3.,4.], [2,2])
        A_false = reshape([1.,2.,3.,4.,5.,6.], [2,3])

        call check(error,is_square(A_true), &
            "is_square(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_square(A_false)), &
            "(.not. is_square(A_false)) failed.")
        if (allocated(error)) return
    end subroutine test_is_square_rsp
    subroutine test_is_square_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true(2,2),A_false(2,3)
        A_true = reshape([1.,2.,3.,4.], [2,2])
        A_false = reshape([1.,2.,3.,4.,5.,6.], [2,3])

        call check(error,is_square(A_true), &
            "is_square(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_square(A_false)), &
            "(.not. is_square(A_false)) failed.")
        if (allocated(error)) return
    end subroutine test_is_square_rdp
    subroutine test_is_square_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true(2,2),A_false(2,3)
        A_true = reshape([1.,2.,3.,4.], [2,2])
        A_false = reshape([1.,2.,3.,4.,5.,6.], [2,3])

        call check(error,is_square(A_true), &
            "is_square(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_square(A_false)), &
            "(.not. is_square(A_false)) failed.")
        if (allocated(error)) return
    end subroutine test_is_square_rqp
    subroutine test_is_square_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true(2,2),A_false(2,3)
        A_true = reshape([cmplx(1.,0.),cmplx(2.,1.),cmplx(3.,0.),cmplx(4.,1.)], [2,2])
        A_false = reshape([cmplx(1.,0.),cmplx(2.,1.),cmplx(3.,0.), &
            cmplx(4.,1.),cmplx(5.,0.),cmplx(6.,1.)], [2,3])

        call check(error,is_square(A_true), &
            "is_square(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_square(A_false)), &
            "(.not. is_square(A_false)) failed.")
        if (allocated(error)) return
    end subroutine test_is_square_csp
    subroutine test_is_square_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true(2,2),A_false(2,3)
        A_true = reshape([cmplx(1.,0.),cmplx(2.,1.),cmplx(3.,0.),cmplx(4.,1.)], [2,2])
        A_false = reshape([cmplx(1.,0.),cmplx(2.,1.),cmplx(3.,0.), &
            cmplx(4.,1.),cmplx(5.,0.),cmplx(6.,1.)], [2,3])

        call check(error,is_square(A_true), &
            "is_square(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_square(A_false)), &
            "(.not. is_square(A_false)) failed.")
        if (allocated(error)) return
    end subroutine test_is_square_cdp
    subroutine test_is_square_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true(2,2),A_false(2,3)
        A_true = reshape([cmplx(1.,0.),cmplx(2.,1.),cmplx(3.,0.),cmplx(4.,1.)], [2,2])
        A_false = reshape([cmplx(1.,0.),cmplx(2.,1.),cmplx(3.,0.), &
            cmplx(4.,1.),cmplx(5.,0.),cmplx(6.,1.)], [2,3])

        call check(error,is_square(A_true), &
            "is_square(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_square(A_false)), &
            "(.not. is_square(A_false)) failed.")
        if (allocated(error)) return
    end subroutine test_is_square_cqp

    !is_diagonal
    subroutine test_is_diagonal_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true_s(2,2),A_false_s(2,2) !square matrices
        real(sp) :: A_true_sf(2,3),A_false_sf(2,3) !short and fat matrices
        real(sp) :: A_true_ts(3,2),A_false_ts(3,2) !tall and skinny matrices
        A_true_s = reshape([1.,0.,0.,4.], [2,2])
        A_false_s = reshape([1.,0.,3.,4.], [2,2])
        A_true_sf = reshape([1.,0.,0.,4.,0.,0.], [2,3])
        A_false_sf = reshape([1.,0.,3.,4.,0.,0.], [2,3])
        A_true_ts = reshape([1.,0.,0.,0.,5.,0.], [3,2])
        A_false_ts = reshape([1.,0.,0.,0.,5.,6.], [3,2])

        call check(error,is_diagonal(A_true_s), &
            "is_diagonal(A_true_s) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_s)), &
            "(.not. is_diagonal(A_false_s)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_sf), &
            "is_diagonal(A_true_sf) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_sf)), &
            "(.not. is_diagonal(A_false_sf)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_ts), &
            "is_diagonal(A_true_ts) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_ts)), &
            "(.not. is_diagonal(A_false_ts)) failed.")
        if (allocated(error)) return
    end subroutine test_is_diagonal_rsp
    subroutine test_is_diagonal_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true_s(2,2),A_false_s(2,2) !square matrices
        real(dp) :: A_true_sf(2,3),A_false_sf(2,3) !short and fat matrices
        real(dp) :: A_true_ts(3,2),A_false_ts(3,2) !tall and skinny matrices
        A_true_s = reshape([1.,0.,0.,4.], [2,2])
        A_false_s = reshape([1.,0.,3.,4.], [2,2])
        A_true_sf = reshape([1.,0.,0.,4.,0.,0.], [2,3])
        A_false_sf = reshape([1.,0.,3.,4.,0.,0.], [2,3])
        A_true_ts = reshape([1.,0.,0.,0.,5.,0.], [3,2])
        A_false_ts = reshape([1.,0.,0.,0.,5.,6.], [3,2])

        call check(error,is_diagonal(A_true_s), &
            "is_diagonal(A_true_s) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_s)), &
            "(.not. is_diagonal(A_false_s)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_sf), &
            "is_diagonal(A_true_sf) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_sf)), &
            "(.not. is_diagonal(A_false_sf)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_ts), &
            "is_diagonal(A_true_ts) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_ts)), &
            "(.not. is_diagonal(A_false_ts)) failed.")
        if (allocated(error)) return
    end subroutine test_is_diagonal_rdp
    subroutine test_is_diagonal_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true_s(2,2),A_false_s(2,2) !square matrices
        real(qp) :: A_true_sf(2,3),A_false_sf(2,3) !short and fat matrices
        real(qp) :: A_true_ts(3,2),A_false_ts(3,2) !tall and skinny matrices
        A_true_s = reshape([1.,0.,0.,4.], [2,2])
        A_false_s = reshape([1.,0.,3.,4.], [2,2])
        A_true_sf = reshape([1.,0.,0.,4.,0.,0.], [2,3])
        A_false_sf = reshape([1.,0.,3.,4.,0.,0.], [2,3])
        A_true_ts = reshape([1.,0.,0.,0.,5.,0.], [3,2])
        A_false_ts = reshape([1.,0.,0.,0.,5.,6.], [3,2])

        call check(error,is_diagonal(A_true_s), &
            "is_diagonal(A_true_s) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_s)), &
            "(.not. is_diagonal(A_false_s)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_sf), &
            "is_diagonal(A_true_sf) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_sf)), &
            "(.not. is_diagonal(A_false_sf)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_ts), &
            "is_diagonal(A_true_ts) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_ts)), &
            "(.not. is_diagonal(A_false_ts)) failed.")
        if (allocated(error)) return
    end subroutine test_is_diagonal_rqp
    subroutine test_is_diagonal_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true_s(2,2),A_false_s(2,2) !square matrices
        complex(sp) :: A_true_sf(2,3),A_false_sf(2,3) !short and fat matrices
        complex(sp) :: A_true_ts(3,2),A_false_ts(3,2) !tall and skinny matrices
        A_true_s = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(4.,1.)], [2,2])
        A_false_s = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,1.)], [2,2])
        A_true_sf = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(4.,1.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_false_sf = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,1.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_true_ts = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(0.,0.)], [3,2])
        A_false_ts = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(6.,1.)], [3,2])

        call check(error,is_diagonal(A_true_s), &
            "is_diagonal(A_true_s) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_s)), &
            "(.not. is_diagonal(A_false_s)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_sf), &
            "is_diagonal(A_true_sf) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_sf)), &
            "(.not. is_diagonal(A_false_sf)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_ts), &
            "is_diagonal(A_true_ts) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_ts)), &
            "(.not. is_diagonal(A_false_ts)) failed.")
        if (allocated(error)) return
    end subroutine test_is_diagonal_csp
    subroutine test_is_diagonal_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true_s(2,2),A_false_s(2,2) !square matrices
        complex(dp) :: A_true_sf(2,3),A_false_sf(2,3) !short and fat matrices
        complex(dp) :: A_true_ts(3,2),A_false_ts(3,2) !tall and skinny matrices
        A_true_s = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(4.,1.)], [2,2])
        A_false_s = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,1.)], [2,2])
        A_true_sf = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(4.,1.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_false_sf = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,1.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_true_ts = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(0.,0.)], [3,2])
        A_false_ts = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(6.,1.)], [3,2])

        call check(error,is_diagonal(A_true_s), &
            "is_diagonal(A_true_s) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_s)), &
            "(.not. is_diagonal(A_false_s)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_sf), &
            "is_diagonal(A_true_sf) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_sf)), &
            "(.not. is_diagonal(A_false_sf)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_ts), &
            "is_diagonal(A_true_ts) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_ts)), &
            "(.not. is_diagonal(A_false_ts)) failed.")
        if (allocated(error)) return
    end subroutine test_is_diagonal_cdp
    subroutine test_is_diagonal_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true_s(2,2),A_false_s(2,2) !square matrices
        complex(qp) :: A_true_sf(2,3),A_false_sf(2,3) !short and fat matrices
        complex(qp) :: A_true_ts(3,2),A_false_ts(3,2) !tall and skinny matrices
        A_true_s = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(4.,1.)], [2,2])
        A_false_s = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,1.)], [2,2])
        A_true_sf = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(4.,1.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_false_sf = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,1.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_true_ts = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(0.,0.)], [3,2])
        A_false_ts = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(6.,1.)], [3,2])

        call check(error,is_diagonal(A_true_s), &
            "is_diagonal(A_true_s) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_s)), &
            "(.not. is_diagonal(A_false_s)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_sf), &
            "is_diagonal(A_true_sf) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_sf)), &
            "(.not. is_diagonal(A_false_sf)) failed.")
        if (allocated(error)) return
        call check(error,is_diagonal(A_true_ts), &
            "is_diagonal(A_true_ts) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_diagonal(A_false_ts)), &
            "(.not. is_diagonal(A_false_ts)) failed.")
        if (allocated(error)) return
    end subroutine test_is_diagonal_cqp

    !is_symmetric
    subroutine test_is_symmetric_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([1.,2.,2.,4.], [2,2])
        A_false_1 = reshape([1.,2.,3.,4.], [2,2])
        A_false_2 = reshape([1.,2.,3.,2.,5.,6.], [3,2]) !nonsquare matrix

        call check(error,is_symmetric(A_true), &
            "is_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_1)), &
            "(.not. is_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_2)), &
            "(.not. is_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_symmetric_rsp
    subroutine test_is_symmetric_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([1.,2.,2.,4.], [2,2])
        A_false_1 = reshape([1.,2.,3.,4.], [2,2])
        A_false_2 = reshape([1.,2.,3.,2.,5.,6.], [3,2]) !nonsquare matrix

        call check(error,is_symmetric(A_true), &
            "is_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_1)), &
            "(.not. is_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_2)), &
            "(.not. is_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_symmetric_rdp
    subroutine test_is_symmetric_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([1.,2.,2.,4.], [2,2])
        A_false_1 = reshape([1.,2.,3.,4.], [2,2])
        A_false_2 = reshape([1.,2.,3.,2.,5.,6.], [3,2]) !nonsquare matrix

        call check(error,is_symmetric(A_true), &
            "is_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_1)), &
            "(.not. is_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_2)), &
            "(.not. is_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_symmetric_rqp
    subroutine test_is_symmetric_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([cmplx(1.,1.),cmplx(2.,1.), &
            cmplx(2.,1.),cmplx(4.,1.)], [2,2])
        A_false_1 = reshape([cmplx(1.,1.),cmplx(2.,1.), &
            cmplx(3.,1.),cmplx(4.,1.)], [2,2])
        A_false_2 = reshape([cmplx(1.,1.),cmplx(2.,1.),cmplx(3.,1.), &
            cmplx(2.,1.),cmplx(5.,1.),cmplx(6.,2.)], [3,2]) !nonsquare matrix

        call check(error,is_symmetric(A_true), &
            "is_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_1)), &
            "(.not. is_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_2)), &
            "(.not. is_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_symmetric_csp
    subroutine test_is_symmetric_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([cmplx(1.,1.),cmplx(2.,1.), &
            cmplx(2.,1.),cmplx(4.,1.)], [2,2])
        A_false_1 = reshape([cmplx(1.,1.),cmplx(2.,1.), &
            cmplx(3.,1.),cmplx(4.,1.)], [2,2])
        A_false_2 = reshape([cmplx(1.,1.),cmplx(2.,1.),cmplx(3.,1.), &
            cmplx(2.,1.),cmplx(5.,1.),cmplx(6.,2.)], [3,2]) !nonsquare matrix

        call check(error,is_symmetric(A_true), &
            "is_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_1)), &
            "(.not. is_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_2)), &
            "(.not. is_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_symmetric_cdp
    subroutine test_is_symmetric_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([cmplx(1.,1.),cmplx(2.,1.), &
            cmplx(2.,1.),cmplx(4.,1.)], [2,2])
        A_false_1 = reshape([cmplx(1.,1.),cmplx(2.,1.), &
            cmplx(3.,1.),cmplx(4.,1.)], [2,2])
        A_false_2 = reshape([cmplx(1.,1.),cmplx(2.,1.),cmplx(3.,1.), &
            cmplx(2.,1.),cmplx(5.,1.),cmplx(6.,2.)], [3,2]) !nonsquare matrix

        call check(error,is_symmetric(A_true), &
            "is_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_1)), &
            "(.not. is_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_symmetric(A_false_2)), &
            "(.not. is_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_symmetric_cqp

    !is_skew_symmetric
    subroutine test_is_skew_symmetric_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2),A_false_diag(2,2)
        A_true = reshape([0.,2.,-2.,0.], [2,2])
        A_false_1 = reshape([0.,2.,-3.,0.], [2,2])
        A_false_2 = reshape([0.,2.,3.,-2.,0.,6.], [3,2]) !nonsquare matrix
        A_false_diag = reshape([1.,2.,-2.,0.], [2,2]) ! non-zero diagonal

        call check(error,is_skew_symmetric(A_true), &
            "is_skew_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_1)), &
            "(.not. is_skew_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_2)), &
            "(.not. is_skew_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_diag)), &
            "(.not. is_skew_symmetric(A_false_diag)) failed.")
        if (allocated(error)) return
    end subroutine test_is_skew_symmetric_rsp
    subroutine test_is_skew_symmetric_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2),A_false_diag(2,2)
        A_true = reshape([0.,2.,-2.,0.], [2,2])
        A_false_1 = reshape([0.,2.,-3.,0.], [2,2])
        A_false_2 = reshape([0.,2.,3.,-2.,0.,6.], [3,2]) !nonsquare matrix
        A_false_diag = reshape([1.,2.,-2.,0.], [2,2]) ! non-zero diagonal

        call check(error,is_skew_symmetric(A_true), &
            "is_skew_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_1)), &
            "(.not. is_skew_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_2)), &
            "(.not. is_skew_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_diag)), &
            "(.not. is_skew_symmetric(A_false_diag)) failed.")
        if (allocated(error)) return
    end subroutine test_is_skew_symmetric_rdp
    subroutine test_is_skew_symmetric_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2),A_false_diag(2,2)
        A_true = reshape([0.,2.,-2.,0.], [2,2])
        A_false_1 = reshape([0.,2.,-3.,0.], [2,2])
        A_false_2 = reshape([0.,2.,3.,-2.,0.,6.], [3,2]) !nonsquare matrix
        A_false_diag = reshape([1.,2.,-2.,0.], [2,2]) ! non-zero diagonal

        call check(error,is_skew_symmetric(A_true), &
            "is_skew_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_1)), &
            "(.not. is_skew_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_2)), &
            "(.not. is_skew_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_diag)), &
            "(.not. is_skew_symmetric(A_false_diag)) failed.")
        if (allocated(error)) return
    end subroutine test_is_skew_symmetric_rqp
    subroutine test_is_skew_symmetric_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2),A_false_diag(2,2)
        A_true = reshape([cmplx(0.,0.),cmplx(2.,1.), &
            -cmplx(2.,1.),cmplx(0.,0.)], [2,2])
        A_false_1 = reshape([cmplx(0.,0.),cmplx(2.,1.), &
            -cmplx(3.,1.),cmplx(0.,0.)], [2,2])
        A_false_2 = reshape([cmplx(0.,0.),cmplx(2.,1.),cmplx(3.,0.), &
            -cmplx(2.,1.),cmplx(0.,0.),cmplx(6.,0.)], [3,2]) !nonsquare matrix
        A_false_diag = reshape([cmplx(0.,1.),cmplx(2.,1.), &
            -cmplx(2.,1.),cmplx(0.,0.)], [2,2]) ! purely imaginary diagonal

        call check(error,is_skew_symmetric(A_true), &
            "is_skew_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_1)), &
            "(.not. is_skew_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_2)), &
            "(.not. is_skew_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_diag)), &
            "(.not. is_skew_symmetric(A_false_diag)) failed.")
        if (allocated(error)) return
    end subroutine test_is_skew_symmetric_csp
    subroutine test_is_skew_symmetric_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2),A_false_diag(2,2)
        A_true = reshape([cmplx(0.,0.),cmplx(2.,1.), &
            -cmplx(2.,1.),cmplx(0.,0.)], [2,2])
        A_false_1 = reshape([cmplx(0.,0.),cmplx(2.,1.), &
            -cmplx(3.,1.),cmplx(0.,0.)], [2,2])
        A_false_2 = reshape([cmplx(0.,0.),cmplx(2.,1.),cmplx(3.,0.), &
            -cmplx(2.,1.),cmplx(0.,0.),cmplx(6.,0.)], [3,2]) !nonsquare matrix
        A_false_diag = reshape([cmplx(0.,1.),cmplx(2.,1.), &
            -cmplx(2.,1.),cmplx(0.,0.)], [2,2]) ! purely imaginary diagonal

        call check(error,is_skew_symmetric(A_true), &
            "is_skew_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_1)), &
            "(.not. is_skew_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_2)), &
            "(.not. is_skew_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_diag)), &
            "(.not. is_skew_symmetric(A_false_diag)) failed.")
        if (allocated(error)) return
    end subroutine test_is_skew_symmetric_cdp
    subroutine test_is_skew_symmetric_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2),A_false_diag(2,2)
        A_true = reshape([cmplx(0.,0.),cmplx(2.,1.), &
            -cmplx(2.,1.),cmplx(0.,0.)], [2,2])
        A_false_1 = reshape([cmplx(0.,0.),cmplx(2.,1.), &
            -cmplx(3.,1.),cmplx(0.,0.)], [2,2])
        A_false_2 = reshape([cmplx(0.,0.),cmplx(2.,1.),cmplx(3.,0.), &
            -cmplx(2.,1.),cmplx(0.,0.),cmplx(6.,0.)], [3,2]) !nonsquare matrix
        A_false_diag = reshape([cmplx(0.,1.),cmplx(2.,1.), &
            -cmplx(2.,1.),cmplx(0.,0.)], [2,2]) ! purely imaginary diagonal

        call check(error,is_skew_symmetric(A_true), &
            "is_skew_symmetric(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_1)), &
            "(.not. is_skew_symmetric(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_2)), &
            "(.not. is_skew_symmetric(A_false_2)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_skew_symmetric(A_false_diag)), &
            "(.not. is_skew_symmetric(A_false_diag)) failed.")
        if (allocated(error)) return
    end subroutine test_is_skew_symmetric_cqp

    !is_hermitian
    subroutine test_is_hermitian_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([1.,2.,2.,4.], [2,2])
        A_false_1 = reshape([1.,2.,3.,4.], [2,2])
        A_false_2 = reshape([1.,2.,3.,2.,5.,6.], [3,2]) !nonsquare matrix

        call check(error,is_hermitian(A_true), &
            "is_hermitian(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_1)), &
            "(.not. is_hermitian(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_2)), &
            "(.not. is_hermitian(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_hermitian_rsp
    subroutine test_is_hermitian_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([1.,2.,2.,4.], [2,2])
        A_false_1 = reshape([1.,2.,3.,4.], [2,2])
        A_false_2 = reshape([1.,2.,3.,2.,5.,6.], [3,2]) !nonsquare matrix

        call check(error,is_hermitian(A_true), &
            "is_hermitian(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_1)), &
            "(.not. is_hermitian(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_2)), &
            "(.not. is_hermitian(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_hermitian_rdp
    subroutine test_is_hermitian_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([1.,2.,2.,4.], [2,2])
        A_false_1 = reshape([1.,2.,3.,4.], [2,2])
        A_false_2 = reshape([1.,2.,3.,2.,5.,6.], [3,2]) !nonsquare matrix

        call check(error,is_hermitian(A_true), &
            "is_hermitian(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_1)), &
            "(.not. is_hermitian(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_2)), &
            "(.not. is_hermitian(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_hermitian_rqp
    subroutine test_is_hermitian_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([cmplx(1.,0.),cmplx(2.,-1.), &
            cmplx(2.,1.),cmplx(4.,0.)], [2,2])
        A_false_1 = reshape([cmplx(1.,0.),cmplx(2.,-1.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_false_2 = reshape([cmplx(1.,0.),cmplx(2.,-1.),cmplx(3.,-1.), &
            cmplx(2.,1.),cmplx(5.,0.),cmplx(6.,-1.)], [3,2]) !nonsquare matrix

        call check(error,is_hermitian(A_true), &
            "is_hermitian(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_1)), &
            "(.not. is_hermitian(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_2)), &
            "(.not. is_hermitian(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_hermitian_csp
    subroutine test_is_hermitian_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([cmplx(1.,0.),cmplx(2.,-1.), &
            cmplx(2.,1.),cmplx(4.,0.)], [2,2])
        A_false_1 = reshape([cmplx(1.,0.),cmplx(2.,-1.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_false_2 = reshape([cmplx(1.,0.),cmplx(2.,-1.),cmplx(3.,-1.), &
            cmplx(2.,1.),cmplx(5.,0.),cmplx(6.,-1.)], [3,2]) !nonsquare matrix

        call check(error,is_hermitian(A_true), &
            "is_hermitian(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_1)), &
            "(.not. is_hermitian(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_2)), &
            "(.not. is_hermitian(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_hermitian_cdp
    subroutine test_is_hermitian_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true(2,2),A_false_1(2,2),A_false_2(3,2)
        A_true = reshape([cmplx(1.,0.),cmplx(2.,-1.), &
            cmplx(2.,1.),cmplx(4.,0.)], [2,2])
        A_false_1 = reshape([cmplx(1.,0.),cmplx(2.,-1.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_false_2 = reshape([cmplx(1.,0.),cmplx(2.,-1.),cmplx(3.,-1.), &
            cmplx(2.,1.),cmplx(5.,0.),cmplx(6.,-1.)], [3,2]) !nonsquare matrix

        call check(error,is_hermitian(A_true), &
            "is_hermitian(A_true) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_1)), &
            "(.not. is_hermitian(A_false_1)) failed.")
        if (allocated(error)) return
        call check(error, (.not. is_hermitian(A_false_2)), &
            "(.not. is_hermitian(A_false_2)) failed.")
        if (allocated(error)) return
    end subroutine test_is_hermitian_cqp

    !is_triangular
    subroutine test_is_triangular_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true_s_u(2,2),A_false_s_u(2,2) !square matrices (upper triangular)
        real(sp) :: A_true_sf_u(2,3),A_false_sf_u(2,3) !short and fat matrices
        real(sp) :: A_true_ts_u(3,2),A_false_ts_u(3,2) !tall and skinny matrices
        real(sp) :: A_true_s_l(2,2),A_false_s_l(2,2) !square matrices (lower triangular)
        real(sp) :: A_true_sf_l(2,3),A_false_sf_l(2,3) !short and fat matrices
        real(sp) :: A_true_ts_l(3,2),A_false_ts_l(3,2) !tall and skinny matrices
        !upper triangular
        A_true_s_u = reshape([1.,0.,3.,4.], [2,2])
        A_false_s_u = reshape([1.,2.,0.,4.], [2,2])
        A_true_sf_u = reshape([1.,0.,3.,4.,0.,6.], [2,3])
        A_false_sf_u = reshape([1.,2.,3.,4.,0.,6.], [2,3])
        A_true_ts_u = reshape([1.,0.,0.,4.,5.,0.], [3,2])
        A_false_ts_u = reshape([1.,0.,0.,4.,5.,6.], [3,2])
        !lower triangular
        A_true_s_l = reshape([1.,2.,0.,4.], [2,2])
        A_false_s_l = reshape([1.,0.,3.,4.], [2,2])
        A_true_sf_l = reshape([1.,2.,0.,4.,0.,0.], [2,3])
        A_false_sf_l = reshape([1.,2.,3.,4.,0.,0.], [2,3])
        A_true_ts_l = reshape([1.,2.,3.,0.,5.,6.], [3,2])
        A_false_ts_l = reshape([1.,2.,3.,4.,5.,6.], [3,2])

        !upper triangular checks
        call check(error,is_triangular(A_true_s_u,'u'), &
            "is_triangular(A_true_s_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_u,'u')), &
            "(.not. is_triangular(A_false_s_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_u,'u'), &
            "is_triangular(A_true_sf_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_u,'u')), &
            "(.not. is_triangular(A_false_sf_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_u,'u'), &
            "is_triangular(A_true_ts_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_u,'u')), &
            "(.not. is_triangular(A_false_ts_u,'u')) failed.")
        if (allocated(error)) return
        !lower triangular checks
        call check(error,is_triangular(A_true_s_l,'l'), &
            "is_triangular(A_true_s_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_l,'l')), &
            "(.not. is_triangular(A_false_s_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_l,'l'), &
            "is_triangular(A_true_sf_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_l,'l')), &
            "(.not. is_triangular(A_false_sf_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_l,'l'), &
            "is_triangular(A_true_ts_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_l,'l')), &
            "(.not. is_triangular(A_false_ts_l,'l')) failed.")
        if (allocated(error)) return
    end subroutine test_is_triangular_rsp
    subroutine test_is_triangular_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true_s_u(2,2),A_false_s_u(2,2) !square matrices (upper triangular)
        real(dp) :: A_true_sf_u(2,3),A_false_sf_u(2,3) !short and fat matrices
        real(dp) :: A_true_ts_u(3,2),A_false_ts_u(3,2) !tall and skinny matrices
        real(dp) :: A_true_s_l(2,2),A_false_s_l(2,2) !square matrices (lower triangular)
        real(dp) :: A_true_sf_l(2,3),A_false_sf_l(2,3) !short and fat matrices
        real(dp) :: A_true_ts_l(3,2),A_false_ts_l(3,2) !tall and skinny matrices
        !upper triangular
        A_true_s_u = reshape([1.,0.,3.,4.], [2,2])
        A_false_s_u = reshape([1.,2.,0.,4.], [2,2])
        A_true_sf_u = reshape([1.,0.,3.,4.,0.,6.], [2,3])
        A_false_sf_u = reshape([1.,2.,3.,4.,0.,6.], [2,3])
        A_true_ts_u = reshape([1.,0.,0.,4.,5.,0.], [3,2])
        A_false_ts_u = reshape([1.,0.,0.,4.,5.,6.], [3,2])
        !lower triangular
        A_true_s_l = reshape([1.,2.,0.,4.], [2,2])
        A_false_s_l = reshape([1.,0.,3.,4.], [2,2])
        A_true_sf_l = reshape([1.,2.,0.,4.,0.,0.], [2,3])
        A_false_sf_l = reshape([1.,2.,3.,4.,0.,0.], [2,3])
        A_true_ts_l = reshape([1.,2.,3.,0.,5.,6.], [3,2])
        A_false_ts_l = reshape([1.,2.,3.,4.,5.,6.], [3,2])

        !upper triangular checks
        call check(error,is_triangular(A_true_s_u,'u'), &
            "is_triangular(A_true_s_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_u,'u')), &
            "(.not. is_triangular(A_false_s_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_u,'u'), &
            "is_triangular(A_true_sf_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_u,'u')), &
            "(.not. is_triangular(A_false_sf_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_u,'u'), &
            "is_triangular(A_true_ts_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_u,'u')), &
            "(.not. is_triangular(A_false_ts_u,'u')) failed.")
        if (allocated(error)) return
        !lower triangular checks
        call check(error,is_triangular(A_true_s_l,'l'), &
            "is_triangular(A_true_s_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_l,'l')), &
            "(.not. is_triangular(A_false_s_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_l,'l'), &
            "is_triangular(A_true_sf_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_l,'l')), &
            "(.not. is_triangular(A_false_sf_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_l,'l'), &
            "is_triangular(A_true_ts_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_l,'l')), &
            "(.not. is_triangular(A_false_ts_l,'l')) failed.")
        if (allocated(error)) return
    end subroutine test_is_triangular_rdp
    subroutine test_is_triangular_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true_s_u(2,2),A_false_s_u(2,2) !square matrices (upper triangular)
        real(qp) :: A_true_sf_u(2,3),A_false_sf_u(2,3) !short and fat matrices
        real(qp) :: A_true_ts_u(3,2),A_false_ts_u(3,2) !tall and skinny matrices
        real(qp) :: A_true_s_l(2,2),A_false_s_l(2,2) !square matrices (lower triangular)
        real(qp) :: A_true_sf_l(2,3),A_false_sf_l(2,3) !short and fat matrices
        real(qp) :: A_true_ts_l(3,2),A_false_ts_l(3,2) !tall and skinny matrices
        !upper triangular
        A_true_s_u = reshape([1.,0.,3.,4.], [2,2])
        A_false_s_u = reshape([1.,2.,0.,4.], [2,2])
        A_true_sf_u = reshape([1.,0.,3.,4.,0.,6.], [2,3])
        A_false_sf_u = reshape([1.,2.,3.,4.,0.,6.], [2,3])
        A_true_ts_u = reshape([1.,0.,0.,4.,5.,0.], [3,2])
        A_false_ts_u = reshape([1.,0.,0.,4.,5.,6.], [3,2])
        !lower triangular
        A_true_s_l = reshape([1.,2.,0.,4.], [2,2])
        A_false_s_l = reshape([1.,0.,3.,4.], [2,2])
        A_true_sf_l = reshape([1.,2.,0.,4.,0.,0.], [2,3])
        A_false_sf_l = reshape([1.,2.,3.,4.,0.,0.], [2,3])
        A_true_ts_l = reshape([1.,2.,3.,0.,5.,6.], [3,2])
        A_false_ts_l = reshape([1.,2.,3.,4.,5.,6.], [3,2])

        !upper triangular checks
        call check(error,is_triangular(A_true_s_u,'u'), &
            "is_triangular(A_true_s_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_u,'u')), &
            "(.not. is_triangular(A_false_s_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_u,'u'), &
            "is_triangular(A_true_sf_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_u,'u')), &
            "(.not. is_triangular(A_false_sf_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_u,'u'), &
            "is_triangular(A_true_ts_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_u,'u')), &
            "(.not. is_triangular(A_false_ts_u,'u')) failed.")
        if (allocated(error)) return
        !lower triangular checks
        call check(error,is_triangular(A_true_s_l,'l'), &
            "is_triangular(A_true_s_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_l,'l')), &
            "(.not. is_triangular(A_false_s_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_l,'l'), &
            "is_triangular(A_true_sf_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_l,'l')), &
            "(.not. is_triangular(A_false_sf_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_l,'l'), &
            "is_triangular(A_true_ts_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_l,'l')), &
            "(.not. is_triangular(A_false_ts_l,'l')) failed.")
        if (allocated(error)) return
    end subroutine test_is_triangular_rqp
    subroutine test_is_triangular_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true_s_u(2,2),A_false_s_u(2,2) !square matrices (upper triangular)
        complex(sp) :: A_true_sf_u(2,3),A_false_sf_u(2,3) !short and fat matrices
        complex(sp) :: A_true_ts_u(3,2),A_false_ts_u(3,2) !tall and skinny matrices
        complex(sp) :: A_true_s_l(2,2),A_false_s_l(2,2) !square matrices (lower triangular)
        complex(sp) :: A_true_sf_l(2,3),A_false_sf_l(2,3) !short and fat matrices
        complex(sp) :: A_true_ts_l(3,2),A_false_ts_l(3,2) !tall and skinny matrices
        !upper triangular
        A_true_s_u = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_false_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.)], [2,2])
        A_true_sf_u = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(6.,0.)], [2,3])
        A_false_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(6.,0.)], [2,3])
        A_true_ts_u = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(0.,0.)], [3,2])
        A_false_ts_u = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])
        !lower triangular
        A_true_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.)], [2,2])
        A_false_s_l = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_true_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_false_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_true_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])
        A_false_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])

        !upper triangular checks
        call check(error,is_triangular(A_true_s_u,'u'), &
            "is_triangular(A_true_s_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_u,'u')), &
            "(.not. is_triangular(A_false_s_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_u,'u'), &
            "is_triangular(A_true_sf_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_u,'u')), &
            "(.not. is_triangular(A_false_sf_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_u,'u'), &
            "is_triangular(A_true_ts_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_u,'u')), &
            "(.not. is_triangular(A_false_ts_u,'u')) failed.")
        if (allocated(error)) return
        !lower triangular checks
        call check(error,is_triangular(A_true_s_l,'l'), &
            "is_triangular(A_true_s_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_l,'l')), &
            "(.not. is_triangular(A_false_s_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_l,'l'), &
            "is_triangular(A_true_sf_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_l,'l')), &
            "(.not. is_triangular(A_false_sf_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_l,'l'), &
            "is_triangular(A_true_ts_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_l,'l')), &
            "(.not. is_triangular(A_false_ts_l,'l')) failed.")
        if (allocated(error)) return
    end subroutine test_is_triangular_csp
    subroutine test_is_triangular_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true_s_u(2,2),A_false_s_u(2,2) !square matrices (upper triangular)
        complex(dp) :: A_true_sf_u(2,3),A_false_sf_u(2,3) !short and fat matrices
        complex(dp) :: A_true_ts_u(3,2),A_false_ts_u(3,2) !tall and skinny matrices
        complex(dp) :: A_true_s_l(2,2),A_false_s_l(2,2) !square matrices (lower triangular)
        complex(dp) :: A_true_sf_l(2,3),A_false_sf_l(2,3) !short and fat matrices
        complex(dp) :: A_true_ts_l(3,2),A_false_ts_l(3,2) !tall and skinny matrices
        !upper triangular
        A_true_s_u = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_false_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.)], [2,2])
        A_true_sf_u = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(6.,0.)], [2,3])
        A_false_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(6.,0.)], [2,3])
        A_true_ts_u = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(0.,0.)], [3,2])
        A_false_ts_u = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])
        !lower triangular
        A_true_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.)], [2,2])
        A_false_s_l = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_true_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_false_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_true_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])
        A_false_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])

        !upper triangular checks
        call check(error,is_triangular(A_true_s_u,'u'), &
            "is_triangular(A_true_s_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_u,'u')), &
            "(.not. is_triangular(A_false_s_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_u,'u'), &
            "is_triangular(A_true_sf_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_u,'u')), &
            "(.not. is_triangular(A_false_sf_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_u,'u'), &
            "is_triangular(A_true_ts_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_u,'u')), &
            "(.not. is_triangular(A_false_ts_u,'u')) failed.")
        if (allocated(error)) return
        !lower triangular checks
        call check(error,is_triangular(A_true_s_l,'l'), &
            "is_triangular(A_true_s_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_l,'l')), &
            "(.not. is_triangular(A_false_s_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_l,'l'), &
            "is_triangular(A_true_sf_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_l,'l')), &
            "(.not. is_triangular(A_false_sf_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_l,'l'), &
            "is_triangular(A_true_ts_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_l,'l')), &
            "(.not. is_triangular(A_false_ts_l,'l')) failed.")
        if (allocated(error)) return
    end subroutine test_is_triangular_cdp
    subroutine test_is_triangular_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true_s_u(2,2),A_false_s_u(2,2) !square matrices (upper triangular)
        complex(qp) :: A_true_sf_u(2,3),A_false_sf_u(2,3) !short and fat matrices
        complex(qp) :: A_true_ts_u(3,2),A_false_ts_u(3,2) !tall and skinny matrices
        complex(qp) :: A_true_s_l(2,2),A_false_s_l(2,2) !square matrices (lower triangular)
        complex(qp) :: A_true_sf_l(2,3),A_false_sf_l(2,3) !short and fat matrices
        complex(qp) :: A_true_ts_l(3,2),A_false_ts_l(3,2) !tall and skinny matrices
        !upper triangular
        A_true_s_u = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_false_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.)], [2,2])
        A_true_sf_u = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(6.,0.)], [2,3])
        A_false_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(6.,0.)], [2,3])
        A_true_ts_u = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(0.,0.)], [3,2])
        A_false_ts_u = reshape([cmplx(1.,1.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])
        !lower triangular
        A_true_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.)], [2,2])
        A_false_s_l = reshape([cmplx(1.,1.),cmplx(0.,0.), &
            cmplx(3.,1.),cmplx(4.,0.)], [2,2])
        A_true_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(0.,0.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_false_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.), &
            cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(0.,0.),cmplx(0.,0.)], [2,3])
        A_true_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(0.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])
        A_false_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.)], [3,2])

        !upper triangular checks
        call check(error,is_triangular(A_true_s_u,'u'), &
            "is_triangular(A_true_s_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_u,'u')), &
            "(.not. is_triangular(A_false_s_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_u,'u'), &
            "is_triangular(A_true_sf_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_u,'u')), &
            "(.not. is_triangular(A_false_sf_u,'u')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_u,'u'), &
            "is_triangular(A_true_ts_u,'u') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_u,'u')), &
            "(.not. is_triangular(A_false_ts_u,'u')) failed.")
        if (allocated(error)) return
        !lower triangular checks
        call check(error,is_triangular(A_true_s_l,'l'), &
            "is_triangular(A_true_s_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_s_l,'l')), &
            "(.not. is_triangular(A_false_s_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_sf_l,'l'), &
            "is_triangular(A_true_sf_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_sf_l,'l')), &
            "(.not. is_triangular(A_false_sf_l,'l')) failed.")
        if (allocated(error)) return
        call check(error,is_triangular(A_true_ts_l,'l'), &
            "is_triangular(A_true_ts_l,'l') failed.")
        if (allocated(error)) return
        call check(error, (.not. is_triangular(A_false_ts_l,'l')), &
            "(.not. is_triangular(A_false_ts_l,'l')) failed.")
        if (allocated(error)) return
    end subroutine test_is_triangular_cqp

    !is_hessenberg
    subroutine test_is_hessenberg_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp) :: A_true_s_u(3,3),A_false_s_u(3,3) !square matrices (upper hessenberg)
        real(sp) :: A_true_sf_u(3,4),A_false_sf_u(3,4) !short and fat matrices
        real(sp) :: A_true_ts_u(4,3),A_false_ts_u(4,3) !tall and skinny matrices
        real(sp) :: A_true_s_l(3,3),A_false_s_l(3,3) !square matrices (lower hessenberg)
        real(sp) :: A_true_sf_l(3,4),A_false_sf_l(3,4) !short and fat matrices
        real(sp) :: A_true_ts_l(4,3),A_false_ts_l(4,3) !tall and skinny matrices
        !upper hessenberg
        A_true_s_u = reshape([1.,2.,0.,4.,5.,6.,7.,8.,9.], [3,3])
        A_false_s_u = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.], [3,3])
        A_true_sf_u = reshape([1.,2.,0.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [3,4])
        A_false_sf_u = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [3,4])
        A_true_ts_u = reshape([1.,2.,0.,0.,5.,6.,7.,0.,9.,10.,11.,12.], [4,3])
        A_false_ts_u = reshape([1.,2.,3.,0.,5.,6.,7.,0.,9.,10.,11.,12.], [4,3])
        !lower hessenberg
        A_true_s_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.], [3,3])
        A_false_s_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.], [3,3])
        A_true_sf_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.,0.,0.,12.], [3,4])
        A_false_sf_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.,0.,11.,12.], [3,4])
        A_true_ts_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,0.,10.,11.,12.], [4,3])
        A_false_ts_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [4,3])

        !upper hessenberg checks
        call check(error,is_hessenberg(A_true_s_u,'u'), &
            "is_hessenberg(A_true_s_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_u,'u')), &
            "(.not. is_hessenberg(A_false_s_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_sf_u,'u'), &
            "is_hessenberg(A_true_sf_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_u,'u')), &
            "(.not. is_hessenberg(A_false_sf_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_ts_u,'u'), &
            "is_hessenberg(A_true_ts_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_u,'u')), &
            "(.not. is_hessenberg(A_false_ts_u,'u')) failed.")
        !lower hessenberg checks
        call check(error,is_hessenberg(A_true_s_l,'l'), &
            "is_hessenberg(A_true_s_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_l,'l')), &
            "(.not. is_hessenberg(A_false_s_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_sf_l,'l'), &
            "is_hessenberg(A_true_sf_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_l,'l')), &
            "(.not. is_hessenberg(A_false_sf_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_ts_l,'l'), &
            "is_hessenberg(A_true_ts_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_l,'l')), &
            "(.not. is_hessenberg(A_false_ts_l,'l')) failed.")
    end subroutine test_is_hessenberg_rsp
    subroutine test_is_hessenberg_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(dp) :: A_true_s_u(3,3),A_false_s_u(3,3) !square matrices (upper hessenberg)
        real(dp) :: A_true_sf_u(3,4),A_false_sf_u(3,4) !short and fat matrices
        real(dp) :: A_true_ts_u(4,3),A_false_ts_u(4,3) !tall and skinny matrices
        real(dp) :: A_true_s_l(3,3),A_false_s_l(3,3) !square matrices (lower hessenberg)
        real(dp) :: A_true_sf_l(3,4),A_false_sf_l(3,4) !short and fat matrices
        real(dp) :: A_true_ts_l(4,3),A_false_ts_l(4,3) !tall and skinny matrices
        !upper hessenberg
        A_true_s_u = reshape([1.,2.,0.,4.,5.,6.,7.,8.,9.], [3,3])
        A_false_s_u = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.], [3,3])
        A_true_sf_u = reshape([1.,2.,0.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [3,4])
        A_false_sf_u = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [3,4])
        A_true_ts_u = reshape([1.,2.,0.,0.,5.,6.,7.,0.,9.,10.,11.,12.], [4,3])
        A_false_ts_u = reshape([1.,2.,3.,0.,5.,6.,7.,0.,9.,10.,11.,12.], [4,3])
        !lower hessenberg
        A_true_s_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.], [3,3])
        A_false_s_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.], [3,3])
        A_true_sf_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.,0.,0.,12.], [3,4])
        A_false_sf_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.,0.,11.,12.], [3,4])
        A_true_ts_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,0.,10.,11.,12.], [4,3])
        A_false_ts_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [4,3])

        !upper hessenberg checks
        call check(error,is_hessenberg(A_true_s_u,'u'), &
            "is_hessenberg(A_true_s_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_u,'u')), &
            "(.not. is_hessenberg(A_false_s_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_sf_u,'u'), &
            "is_hessenberg(A_true_sf_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_u,'u')), &
            "(.not. is_hessenberg(A_false_sf_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_ts_u,'u'), &
            "is_hessenberg(A_true_ts_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_u,'u')), &
            "(.not. is_hessenberg(A_false_ts_u,'u')) failed.")
        !lower hessenberg checks
        call check(error,is_hessenberg(A_true_s_l,'l'), &
            "is_hessenberg(A_true_s_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_l,'l')), &
            "(.not. is_hessenberg(A_false_s_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_sf_l,'l'), &
            "is_hessenberg(A_true_sf_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_l,'l')), &
            "(.not. is_hessenberg(A_false_sf_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_ts_l,'l'), &
            "is_hessenberg(A_true_ts_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_l,'l')), &
            "(.not. is_hessenberg(A_false_ts_l,'l')) failed.")
    end subroutine test_is_hessenberg_rdp
    subroutine test_is_hessenberg_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(qp) :: A_true_s_u(3,3),A_false_s_u(3,3) !square matrices (upper hessenberg)
        real(qp) :: A_true_sf_u(3,4),A_false_sf_u(3,4) !short and fat matrices
        real(qp) :: A_true_ts_u(4,3),A_false_ts_u(4,3) !tall and skinny matrices
        real(qp) :: A_true_s_l(3,3),A_false_s_l(3,3) !square matrices (lower hessenberg)
        real(qp) :: A_true_sf_l(3,4),A_false_sf_l(3,4) !short and fat matrices
        real(qp) :: A_true_ts_l(4,3),A_false_ts_l(4,3) !tall and skinny matrices
        !upper hessenberg
        A_true_s_u = reshape([1.,2.,0.,4.,5.,6.,7.,8.,9.], [3,3])
        A_false_s_u = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.], [3,3])
        A_true_sf_u = reshape([1.,2.,0.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [3,4])
        A_false_sf_u = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [3,4])
        A_true_ts_u = reshape([1.,2.,0.,0.,5.,6.,7.,0.,9.,10.,11.,12.], [4,3])
        A_false_ts_u = reshape([1.,2.,3.,0.,5.,6.,7.,0.,9.,10.,11.,12.], [4,3])
        !lower hessenberg
        A_true_s_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.], [3,3])
        A_false_s_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.], [3,3])
        A_true_sf_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.,0.,0.,12.], [3,4])
        A_false_sf_l = reshape([1.,2.,3.,4.,5.,6.,0.,8.,9.,0.,11.,12.], [3,4])
        A_true_ts_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,0.,10.,11.,12.], [4,3])
        A_false_ts_l = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.], [4,3])

        !upper hessenberg checks
        call check(error,is_hessenberg(A_true_s_u,'u'), &
            "is_hessenberg(A_true_s_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_u,'u')), &
            "(.not. is_hessenberg(A_false_s_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_sf_u,'u'), &
            "is_hessenberg(A_true_sf_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_u,'u')), &
            "(.not. is_hessenberg(A_false_sf_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_ts_u,'u'), &
            "is_hessenberg(A_true_ts_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_u,'u')), &
            "(.not. is_hessenberg(A_false_ts_u,'u')) failed.")
        !lower hessenberg checks
        call check(error,is_hessenberg(A_true_s_l,'l'), &
            "is_hessenberg(A_true_s_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_l,'l')), &
            "(.not. is_hessenberg(A_false_s_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_sf_l,'l'), &
            "is_hessenberg(A_true_sf_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_l,'l')), &
            "(.not. is_hessenberg(A_false_sf_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_ts_l,'l'), &
            "is_hessenberg(A_true_ts_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_l,'l')), &
            "(.not. is_hessenberg(A_false_ts_l,'l')) failed.")
    end subroutine test_is_hessenberg_rqp
    subroutine test_is_hessenberg_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(sp) :: A_true_s_u(3,3),A_false_s_u(3,3) !square matrices (upper hessenberg)
        complex(sp) :: A_true_sf_u(3,4),A_false_sf_u(3,4) !short and fat matrices
        complex(sp) :: A_true_ts_u(4,3),A_false_ts_u(4,3) !tall and skinny matrices
        complex(sp) :: A_true_s_l(3,3),A_false_s_l(3,3) !square matrices (lower hessenberg)
        complex(sp) :: A_true_sf_l(3,4),A_false_sf_l(3,4) !short and fat matrices
        complex(sp) :: A_true_ts_l(4,3),A_false_ts_l(4,3) !tall and skinny matrices
        !upper hessenberg
        A_true_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_false_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_true_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_false_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_true_ts_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(0.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        A_false_ts_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(0.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(0.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        !lower hessenberg
        A_true_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_false_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_true_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(0.,0.),cmplx(0.,0.),cmplx(12.,0.)], [3,4])
        A_false_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(0.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_true_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(8.,0.), &
            cmplx(0.,0.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        A_false_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(8.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])

        !upper hessenberg checks
        call check(error,is_hessenberg(A_true_s_u,'u'), &
            "is_hessenberg(A_true_s_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_u,'u')), &
            "(.not. is_hessenberg(A_false_s_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_sf_u,'u'), &
            "is_hessenberg(A_true_sf_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_u,'u')), &
            "(.not. is_hessenberg(A_false_sf_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_ts_u,'u'), &
            "is_hessenberg(A_true_ts_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_u,'u')), &
            "(.not. is_hessenberg(A_false_ts_u,'u')) failed.")
        !lower hessenberg checks
        call check(error,is_hessenberg(A_true_s_l,'l'), &
            "is_hessenberg(A_true_s_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_l,'l')), &
            "(.not. is_hessenberg(A_false_s_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_sf_l,'l'), &
            "is_hessenberg(A_true_sf_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_l,'l')), &
            "(.not. is_hessenberg(A_false_sf_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_ts_l,'l'), &
            "is_hessenberg(A_true_ts_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_l,'l')), &
            "(.not. is_hessenberg(A_false_ts_l,'l')) failed.")
    end subroutine test_is_hessenberg_csp
    subroutine test_is_hessenberg_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(dp) :: A_true_s_u(3,3),A_false_s_u(3,3) !square matrices (upper hessenberg)
        complex(dp) :: A_true_sf_u(3,4),A_false_sf_u(3,4) !short and fat matrices
        complex(dp) :: A_true_ts_u(4,3),A_false_ts_u(4,3) !tall and skinny matrices
        complex(dp) :: A_true_s_l(3,3),A_false_s_l(3,3) !square matrices (lower hessenberg)
        complex(dp) :: A_true_sf_l(3,4),A_false_sf_l(3,4) !short and fat matrices
        complex(dp) :: A_true_ts_l(4,3),A_false_ts_l(4,3) !tall and skinny matrices
        !upper hessenberg
        A_true_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_false_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_true_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_false_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_true_ts_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(0.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        A_false_ts_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(0.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(0.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        !lower hessenberg
        A_true_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_false_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_true_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(0.,0.),cmplx(0.,0.),cmplx(12.,0.)], [3,4])
        A_false_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(0.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_true_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(8.,0.), &
            cmplx(0.,0.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        A_false_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(8.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])

        !upper hessenberg checks
        call check(error,is_hessenberg(A_true_s_u,'u'), &
            "is_hessenberg(A_true_s_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_u,'u')), &
            "(.not. is_hessenberg(A_false_s_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_sf_u,'u'), &
            "is_hessenberg(A_true_sf_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_u,'u')), &
            "(.not. is_hessenberg(A_false_sf_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_ts_u,'u'), &
            "is_hessenberg(A_true_ts_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_u,'u')), &
            "(.not. is_hessenberg(A_false_ts_u,'u')) failed.")
        !lower hessenberg checks
        call check(error,is_hessenberg(A_true_s_l,'l'), &
            "is_hessenberg(A_true_s_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_l,'l')), &
            "(.not. is_hessenberg(A_false_s_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_sf_l,'l'), &
            "is_hessenberg(A_true_sf_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_l,'l')), &
            "(.not. is_hessenberg(A_false_sf_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_ts_l,'l'), &
            "is_hessenberg(A_true_ts_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_l,'l')), &
            "(.not. is_hessenberg(A_false_ts_l,'l')) failed.")
    end subroutine test_is_hessenberg_cdp
    subroutine test_is_hessenberg_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        complex(qp) :: A_true_s_u(3,3),A_false_s_u(3,3) !square matrices (upper hessenberg)
        complex(qp) :: A_true_sf_u(3,4),A_false_sf_u(3,4) !short and fat matrices
        complex(qp) :: A_true_ts_u(4,3),A_false_ts_u(4,3) !tall and skinny matrices
        complex(qp) :: A_true_s_l(3,3),A_false_s_l(3,3) !square matrices (lower hessenberg)
        complex(qp) :: A_true_sf_l(3,4),A_false_sf_l(3,4) !short and fat matrices
        complex(qp) :: A_true_ts_l(4,3),A_false_ts_l(4,3) !tall and skinny matrices
        !upper hessenberg
        A_true_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_false_s_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_true_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_false_sf_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_true_ts_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(0.,0.),cmplx(0.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(0.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        A_false_ts_u = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(0.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(0.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        !lower hessenberg
        A_true_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_false_s_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(7.,1.),cmplx(8.,0.),cmplx(9.,1.)], [3,3])
        A_true_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(0.,0.),cmplx(0.,0.),cmplx(12.,0.)], [3,4])
        A_false_sf_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.), &
            cmplx(4.,0.),cmplx(5.,1.),cmplx(6.,0.), &
            cmplx(0.,0.),cmplx(8.,0.),cmplx(9.,1.), &
            cmplx(0.,0.),cmplx(11.,1.),cmplx(12.,0.)], [3,4])
        A_true_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(8.,0.), &
            cmplx(0.,0.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])
        A_false_ts_l = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(3.,1.),cmplx(4.,0.), &
            cmplx(5.,1.),cmplx(6.,0.),cmplx(7.,1.),cmplx(8.,0.), &
            cmplx(9.,1.),cmplx(10.,0.),cmplx(11.,1.),cmplx(12.,0.)], [4,3])

        !upper hessenberg checks
        call check(error,is_hessenberg(A_true_s_u,'u'), &
            "is_hessenberg(A_true_s_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_u,'u')), &
            "(.not. is_hessenberg(A_false_s_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_sf_u,'u'), &
            "is_hessenberg(A_true_sf_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_u,'u')), &
            "(.not. is_hessenberg(A_false_sf_u,'u')) failed.")
        call check(error,is_hessenberg(A_true_ts_u,'u'), &
            "is_hessenberg(A_true_ts_u,'u') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_u,'u')), &
            "(.not. is_hessenberg(A_false_ts_u,'u')) failed.")
        !lower hessenberg checks
        call check(error,is_hessenberg(A_true_s_l,'l'), &
            "is_hessenberg(A_true_s_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_s_l,'l')), &
            "(.not. is_hessenberg(A_false_s_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_sf_l,'l'), &
            "is_hessenberg(A_true_sf_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_sf_l,'l')), &
            "(.not. is_hessenberg(A_false_sf_l,'l')) failed.")
        call check(error,is_hessenberg(A_true_ts_l,'l'), &
            "is_hessenberg(A_true_ts_l,'l') failed.")
        call check(error, (.not. is_hessenberg(A_false_ts_l,'l')), &
            "(.not. is_hessenberg(A_false_ts_l,'l')) failed.")
    end subroutine test_is_hessenberg_cqp

end module

program tester
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg_matrix_property_checks,only:collect_linalg_matrix_property_checks
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_matrix_property_checks",collect_linalg_matrix_property_checks) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program
