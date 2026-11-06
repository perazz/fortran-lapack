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

        call add_test(tests,new_unittest("gesv_info",test_gesv_info))
        call add_test(tests,new_unittest("getrf_info",test_getrf_info))
        call add_test(tests,new_unittest("getri_info",test_getri_info))
        call add_test(tests,new_unittest("gecon_info",test_gecon_info))
        call add_test(tests,new_unittest("posv_info",test_posv_info))
        call add_test(tests,new_unittest("potrs_info",test_potrs_info))
        call add_test(tests,new_unittest("potrf_info",test_potrf_info))
        call add_test(tests,new_unittest("geqrf_info",test_geqrf_info))
        call add_test(tests,new_unittest("geqp3_info",test_geqp3_info))
        call add_test(tests,new_unittest("orgqr_info",test_orgqr_info))
        call add_test(tests,new_unittest("gees_info",test_gees_info))
        call add_test(tests,new_unittest("geev_info",test_geev_info))
        call add_test(tests,new_unittest("ggev_info",test_ggev_info))
        call add_test(tests,new_unittest("heev_info",test_heev_info))
        call add_test(tests,new_unittest("gesdd_info",test_gesdd_info))
        call add_test(tests,new_unittest("gelsd_info",test_gelsd_info))
        call add_test(tests,new_unittest("gglse_info",test_gglse_info))
        call add_test(tests,new_unittest("ggglm_info",test_ggglm_info))

    end subroutine test_info_handlers

    subroutine test_gesv_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-7,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gesv_info(caller,info,nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'gesv',info,err,flag)

    end subroutine test_gesv_info

    subroutine test_getrf_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_getrf_info(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'getrf',info,err,flag)

    end subroutine test_getrf_info

    subroutine test_getri_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-3,-6,1,-2]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_getri_info(caller,info,nmat,nmat,err)
        call check_states(error,'getri',info,err,flag)

    end subroutine test_getri_info

    subroutine test_gecon_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-5,1,-3,2]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,INTERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gecon_info(caller,info,nmat,nmat,err)
        call check_states(error,'gecon',info,err,flag)

    end subroutine test_gecon_info

    subroutine test_posv_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_posv_info(caller,info,'L',nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'posv',info,err,flag)

    end subroutine test_posv_info

    subroutine test_potrs_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,VALERR,VALERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_potrs_info(caller,info,'L',nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'potrs',info,err,flag)

    end subroutine test_potrs_info

    subroutine test_potrf_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_potrf_info(caller,info,'L',nmat,nmat,err)
        call check_states(error,'potrf',info,err,flag)

    end subroutine test_potrf_info

    subroutine test_geqrf_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-7,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,INTERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_geqrf_info(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'geqrf',info,err,flag)

    end subroutine test_geqrf_info

    subroutine test_geqp3_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-8,1,-3]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,INTERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_geqp3_info(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'geqp3',info,err,flag)

    end subroutine test_geqp3_info

    subroutine test_orgqr_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-8,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,INTERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_orgqr_info(caller,info,nmat,nmat,nmat,nmat,err)
        call check_states(error,'orgqr',info,err,flag)

    end subroutine test_orgqr_info

    !> Both argument lists: ldvs is -11 (real) or -10 (complex), lwork -13 or -12
    subroutine test_gees_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-4,-6,-10,-11,-12,-13,1,nmat,nmat+1,nmat+2,nmat+3,-3]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,VALERR,VALERR,INTERR,INTERR, &
                                              ALGERR,ALGERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gees_info(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'gees',info,err,flag)

    end subroutine test_gees_info

    !> Both argument lists: ldvl is -9 (real) or -8 (complex), ldvr -11 or -10, lwork -13 or -12
    subroutine test_geev_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-8,-9,-10,-11,-12,-13,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR, &
                                              INTERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_geev_info(caller,info,nmat,nmat,err)
        call check_states(error,'geev',info,err,flag)

    end subroutine test_geev_info

    !> Both argument lists: ldvl is -12 (real) or -11 (complex), ldvr -14 or -13, lwork -16 or -15
    subroutine test_ggev_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-11,-12,-13,-14,-15,-16,1,nmat+1,nmat+2,nmat+3,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR, &
                                              INTERR,INTERR,ALGERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_ggev_info(caller,info,nmat,nmat,nmat,err)
        call check_states(error,'ggev',info,err,flag)

    end subroutine test_ggev_info

    subroutine test_heev_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-8,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,INTERR,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_heev_info(caller,info,nmat,nmat,err)
        call check_states(error,'heev',info,err,flag)

    end subroutine test_heev_info

    subroutine test_gesdd_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-4,-5,-8,-10,-12,1,-6]
        integer(ilp), parameter :: flag(*) = [OK,INTERR,VALERR,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gesdd_info(caller,info,nmat,nmat,err)
        call check_states(error,'gesdd',info,err,flag)

    end subroutine test_gesdd_info

    subroutine test_gelsd_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-12,1,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gelsd_info(caller,info,nmat,1_ilp,nmat,nmat,err)
        call check_states(error,'gelsd',info,err,flag)

    end subroutine test_gelsd_info

    subroutine test_gglse_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-12,1,2,3,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_gglse_info(caller,info,nmat,nmat,2_ilp,err)
        call check_states(error,'gglse',info,err,flag)

    end subroutine test_gglse_info

    subroutine test_ggglm_info(error)
        type(error_type), allocatable, intent(out) :: error

        integer(ilp), parameter :: info(*) = [0,-1,-2,-3,-5,-7,-12,1,2,3,-4]
        integer(ilp), parameter :: flag(*) = [OK,VALERR,VALERR,VALERR,VALERR,VALERR,INTERR,ALGERR,ALGERR,INTERR,INTERR]
        type(la_state) :: err(size(info))

        call handle_ggglm_info(caller,info,nmat,nmat,err)
        call check_states(error,'ggglm',info,err,flag)

    end subroutine test_ggglm_info

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
