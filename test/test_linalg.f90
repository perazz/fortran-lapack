
module test_linalg
    use testdrive,only:new_unittest,unittest_type,error_type,check
    use la_constants,only:sp,dp,xdp,qp,ilp,lk
    use linear_algebra,only:diag,eye

    implicit none

    real(sp),parameter :: sptol = 1000*epsilon(1._sp)
    real(dp),parameter :: dptol = 1000*epsilon(1._dp)
#ifdef LA_WITH_XDP
    real(xdp),parameter :: xdptol = 1000*epsilon(1._xdp)
#endif
#ifdef LA_WITH_QP
    real(qp),parameter :: qptol = 1000*epsilon(1._qp)
#endif

contains

    !> Collect all exported unit tests
    subroutine collect_linalg(testsuite)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("diag_rsp",test_diag_rsp), &
            new_unittest("diag_rdp",test_diag_rdp), &
#ifdef LA_WITH_XDP
            new_unittest("diag_rxdp",test_diag_rxdp), &
#endif
#ifdef LA_WITH_QP
            new_unittest("diag_rqp",test_diag_rqp), &
#endif
            new_unittest("eye",test_eye) &
            ]

    end subroutine collect_linalg

    subroutine test_eye(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp),allocatable :: rye(:,:)
        complex(sp) :: cye(7,7)
        integer :: i

        call check(error,all(eye(3,3) == diag([(1.0_dp,i=1,3)])), &
            "all(eye(3,3) == diag([(1.0_dp,i=1,3)])) failed.")
        if (allocated(error)) return

        rye = real(eye(3,4),sp)
        call check(error,sum(abs(rye(:,1:3) - diag([(1.0_sp,i=1,3)]))) < sptol, &
            "sum(abs(rye(:,1:3) - diag([(1.0_sp,i=1,3)]))) < sptol failed")
        if (allocated(error)) return

        call check(error,all(eye(5) == diag([(1.0_dp,i=1,5)])), &
            "all(eye(5) == diag([(1.0_dp,i=1,5)] failed.")
        if (allocated(error)) return

        rye = real(eye(6),sp)
        call check(error,sum(rye - diag([(1.0_sp,i=1,6)])) < sptol, &
            "sum(rye - diag([(1.0_sp,i=1,6)])) < sptol failed.")
        if (allocated(error)) return

        cye = real(eye(7),sp)
        call check(error,abs(sum([(cye(i,i),i=1,7)]) - cmplx(7.0_sp,0.0_sp,kind=sp)) < sptol, &
            "abs(sum(diagonal(cye)) - cmplx(7.0_sp,0.0_sp,kind=sp)) < sptol failed.")

    end subroutine test_eye

    subroutine test_diag_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")

    end subroutine test_diag_rsp

    subroutine test_diag_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")

    end subroutine test_diag_rdp

#ifdef LA_WITH_XDP
    subroutine test_diag_rxdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(xdp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")

    end subroutine test_diag_rxdp
#endif

#ifdef LA_WITH_QP
    subroutine test_diag_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")

    end subroutine test_diag_rqp
#endif

end module test_linalg

program test_eye_diag
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg,only:collect_linalg
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg",collect_linalg) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program test_eye_diag
