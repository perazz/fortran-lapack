! Test the shared LAPACK info handlers: every documented code, a positive code and an undocumented one
module test_linalg_info
    use testdrive, only: error_type, check, new_unittest, unittest_type
    use la_constants, only: ilp
    use la_state_type, only: la_state,LINALG_SUCCESS,LINALG_VALUE_ERROR,LINALG_ERROR,LINALG_INTERNAL_ERROR
    use la_info

    implicit none (type,external)
    private

    public :: test_info_handlers

    !> The caller name each handler receives, and must return as where_at
    character(*), parameter :: caller = 'test_caller'

    !> Short names for the expected flags
    integer(ilp), parameter :: OK     = LINALG_SUCCESS
    integer(ilp), parameter :: VALERR = LINALG_VALUE_ERROR
    integer(ilp), parameter :: ALGERR = LINALG_ERROR
    integer(ilp), parameter :: INTERR = LINALG_INTERNAL_ERROR

    !> Matrix size the handlers receive, where the positive codes depend on it
    integer(ilp), parameter :: nmat = 4

    contains

    !> LAPACK info handler tests
    subroutine test_info_handlers(tests)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: tests(:)

        allocate(tests(0))

        call add_test(tests,new_unittest("handle_gesv",test_handle_gesv))
        call add_test(tests,new_unittest("handle_getrf",test_handle_getrf))
        call add_test(tests,new_unittest("handle_getri",test_handle_getri))
        call add_test(tests,new_unittest("handle_gecon",test_handle_gecon))
        call add_test(tests,new_unittest("handle_posv",test_handle_posv))
        call add_test(tests,new_unittest("handle_potrs",test_handle_potrs))
        call add_test(tests,new_unittest("handle_potrf",test_handle_potrf))
        call add_test(tests,new_unittest("handle_geqrf",test_handle_geqrf))
        call add_test(tests,new_unittest("handle_geqp3",test_handle_geqp3))
        call add_test(tests,new_unittest("handle_orgqr",test_handle_orgqr))
        call add_test(tests,new_unittest("handle_gees",test_handle_gees))
        call add_test(tests,new_unittest("handle_geev",test_handle_geev))
        call add_test(tests,new_unittest("handle_ggev",test_handle_ggev))
        call add_test(tests,new_unittest("handle_heev",test_handle_heev))
        call add_test(tests,new_unittest("handle_gesdd",test_handle_gesdd))
        call add_test(tests,new_unittest("handle_gelsd",test_handle_gelsd))
        call add_test(tests,new_unittest("handle_gglse",test_handle_gglse))
        call add_test(tests,new_unittest("handle_ggglm",test_handle_ggglm))

    end subroutine test_info_handlers

    subroutine test_handle_gesv(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-7,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gesv(caller,info,nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'gesv',info,err,flag)

    end subroutine test_handle_gesv

    subroutine test_handle_getrf(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_getrf(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'getrf',info,err,flag)

    end subroutine test_handle_getrf

    subroutine test_handle_getri(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-3,-6,1,-2]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_getri(caller,info,nmat,nmat,err)
        call check_states(error,'getri',info,err,flag)

    end subroutine test_handle_getri

    subroutine test_handle_gecon(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-5,1,-3,2]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,INTERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gecon(caller,info,nmat,nmat,err)
        call check_states(error,'gecon',info,err,flag)

    end subroutine test_handle_gecon

    subroutine test_handle_posv(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_posv(caller,info,'L',nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'posv',info,err,flag)

    end subroutine test_handle_posv

    subroutine test_handle_potrs(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,VALERR,VALERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_potrs(caller,info,'L',nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'potrs',info,err,flag)

    end subroutine test_handle_potrs

    subroutine test_handle_potrf(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_potrf(caller,info,'L',nmat,nmat,err)
        call check_states(error,'potrf',info,err,flag)

    end subroutine test_handle_potrf

    subroutine test_handle_geqrf(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-7,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,INTERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_geqrf(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'geqrf',info,err,flag)

    end subroutine test_handle_geqrf

    subroutine test_handle_geqp3(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-8,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,INTERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_geqp3(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'geqp3',info,err,flag)

    end subroutine test_handle_geqp3

    subroutine test_handle_orgqr(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-8,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,INTERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_orgqr(caller,info,nmat,nmat,nmat,nmat,err)
        call check_states(error,'orgqr',info,err,flag)

    end subroutine test_handle_orgqr

    !> Both argument lists: ldvs is -11 (real) or -10 (complex), lwork -13 or -12
    subroutine test_handle_gees(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-6,-10,-11,-12,-13,1,nmat,nmat+1,nmat+2,nmat+3,-3]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,VALERR,VALERR,INTERR,INTERR, &
                                              ALGERR,ALGERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gees(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'gees',info,err,flag)

    end subroutine test_handle_gees

    !> Both argument lists: ldvl is -9 (real) or -8 (complex), ldvr -11 or -10, lwork -13 or -12
    subroutine test_handle_geev(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-8,-9,-10,-11,-12,-13,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR, &
                                              INTERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_geev(caller,info,nmat,nmat,err)
        call check_states(error,'geev',info,err,flag)

    end subroutine test_handle_geev

    !> Both argument lists: ldvl is -12 (real) or -11 (complex), ldvr -14 or -13, lwork -16 or -15
    subroutine test_handle_ggev(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-11,-12,-13,-14,-15,-16,1,nmat+1,nmat+2,nmat+3,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR, &
                                              INTERR,INTERR,ALGERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_ggev(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'ggev',info,err,flag)

    end subroutine test_handle_ggev

    subroutine test_handle_heev(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-8,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_heev(caller,info,nmat,nmat,err)
        call check_states(error,'heev',info,err,flag)

    end subroutine test_handle_heev

    subroutine test_handle_gesdd(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-4,-5,-8,-10,-12,1,-6]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gesdd(caller,info,nmat,nmat,err)
        call check_states(error,'gesdd',info,err,flag)

    end subroutine test_handle_gesdd

    subroutine test_handle_gelsd(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-12,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gelsd(caller,info,nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'gelsd',info,err,flag)

    end subroutine test_handle_gelsd

    subroutine test_handle_gglse(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-12,1,2,3,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gglse(caller,info,nmat,nmat,2_ilp,err)
        call check_states(error,'gglse',info,err,flag)

    end subroutine test_handle_gglse

    subroutine test_handle_ggglm(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-12,1,2,3,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_ggglm(caller,info,nmat,nmat,err)
        call check_states(error,'ggglm',info,err,flag)

    end subroutine test_handle_ggglm

    !> Each info code returns its flag, and names the caller exactly when it is an error
    subroutine check_states(error,routine,info,err,flag)
        type(error_type), allocatable, intent(out) :: error
        character(*), intent(in) :: routine
        integer(ilp), intent(in) :: info(:)
        type(la_state), intent(in) :: err(:)
        integer(ilp), intent(in) :: flag(:)

        integer :: i
        character(len=64) :: label

        do i = 1, size(info)
            write(label,'(a,1x,a,i0)') routine,'info=',info(i)

            call check(error, err(i)%state, flag(i), trim(label)//' returned '//err(i)%print())
            if (allocated(error)) return

            if (flag(i)==OK) then
                call check(error, len_trim(err(i)%where_at)==0, trim(label)//' set where_at on success')
            else
                call check(error, err(i)%where_at==caller, trim(label)//' returned where_at='//trim(err(i)%where_at))
            end if
            if (allocated(error)) return
        end do

    end subroutine check_states

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type), allocatable, intent(inout) :: tests(:)
        type(unittest_type), intent(in) :: new_test

        integer :: n
        type(unittest_type), allocatable :: new_tests(:)

        if (allocated(tests)) then
            n = size(tests)
        else
            n = 0
        end if

        allocate(new_tests(n+1))
        if (n>0) new_tests(1:n) = tests(1:n)
                 new_tests(1+n) = new_test
        call move_alloc(from=new_tests,to=tests)

    end subroutine add_test

end module test_linalg_info

program test_info
     use, intrinsic :: iso_fortran_env, only : error_unit
     use testdrive, only : run_testsuite, new_testsuite, testsuite_type
     use test_linalg_info, only : test_info_handlers
     implicit none
     integer :: stat, is
     type(testsuite_type), allocatable :: testsuites(:)
     character(len=*), parameter :: fmt = '("#", *(1x, a))'

     stat = 0

     testsuites = [ &
         new_testsuite("linalg_info", test_info_handlers) &
         ]

     do is = 1, size(testsuites)
         write(error_unit, fmt) "Testing:", testsuites(is)%name
         call run_testsuite(testsuites(is)%collect, error_unit, stat)
     end do

     if (stat > 0) then
         write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
         error stop
     end if
end program test_info
