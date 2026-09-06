! Test the precision kinds the build carries
module test_linalg_kinds
    use testdrive,only:error_type,check,new_unittest,unittest_type
    use la_constants

    implicit none(type,external)
    private

    public :: test_kind_parameters

    contains

    !> Precision kind tests
    subroutine test_kind_parameters(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("optional_kinds",test_optional_kinds))

    end subroutine test_kind_parameters

    !> Each optional kind is a usable kind number exactly when its flag says so
    subroutine test_optional_kinds(error)
        type(error_type),allocatable,intent(out) :: error

        call check(error,la_with_qp .eqv. (qp > 0),'la_with_qp disagrees with the qp kind number')
        if (allocated(error)) return

        call check(error,la_with_xdp .eqv. (xdp > 0),'la_with_xdp disagrees with the xdp kind number')
        if (allocated(error)) return

#ifdef LA_WITH_QP
        call check(error,precision(1.0_qp) >= 30,'qp holds fewer than 30 decimal digits')
        if (allocated(error)) return
#endif

#ifdef LA_WITH_XDP
        call check(error,digits(1.0_xdp) > digits(1.0_dp),'xdp is no wider than dp')
        if (allocated(error)) return

        call check(error,radix(1.0_xdp) == 2,'xdp is not a binary kind')
        if (allocated(error)) return
#endif

#if defined(LA_WITH_XDP) && defined(LA_WITH_QP)
        call check(error,digits(1.0_xdp) < digits(1.0_qp),'xdp is no narrower than qp')
        if (allocated(error)) return
#endif

    end subroutine test_optional_kinds

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

end module test_linalg_kinds

program test_kinds
     use,intrinsic :: iso_fortran_env,only:error_unit
     use testdrive,only:run_testsuite,new_testsuite,testsuite_type
     use test_linalg_kinds,only:test_kind_parameters
     implicit none
     integer :: stat,is
     type(testsuite_type),allocatable :: testsuites(:)
     character(len=*),parameter :: fmt = '("#", *(1x, a))'

     stat = 0

     testsuites = [ &
         new_testsuite("linalg_kinds",test_kind_parameters) &
         ]

     do is = 1,size(testsuites)
         write (error_unit,fmt) "Testing:",testsuites(is)%name
         call run_testsuite(testsuites(is)%collect,error_unit,stat)
     end do

     if (stat > 0) then
         write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
         error stop
     end if
end program test_kinds
