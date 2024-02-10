
module test_linalg_solve
    use stdlib_linalg_constants
    use stdlib_linalg_state
    use stdlib_linalg_solve

    implicit none (type,external)

    contains

    !> Solve real linear system
    subroutine test_solve(error)
        logical, intent(out) :: error

        call test_ssolve(error)
        if (error) return
        call test_ssolve_multiple(error)
        if (error) return
        call test_dsolve(error)
        if (error) return
        call test_dsolve_multiple(error)
        if (error) return
        call test_qsolve(error)
        if (error) return
        call test_qsolve_multiple(error)
        if (error) return
        call test_csolve(error)
        if (error) return
        call test_2x2_csolve(error)
        if (error) return
        call test_zsolve(error)
        if (error) return
        call test_2x2_zsolve(error)
        if (error) return
        call test_wsolve(error)
        if (error) return
        call test_2x2_wsolve(error)
        if (error) return
    end subroutine test_solve

    !> Simple linear system
    subroutine test_ssolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(sp) :: A(3,3) = transpose(reshape([real(sp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(sp) :: b  (3) = [real(sp) :: 1, 4, -1]
        real(sp) :: res(3) = [real(sp) :: -2, -2, 3]
        real(sp) :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_sp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_sp)
        print *, 'state = ',state%print()


    end subroutine test_ssolve
    subroutine test_dsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(dp) :: A(3,3) = transpose(reshape([real(dp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(dp) :: b  (3) = [real(dp) :: 1, 4, -1]
        real(dp) :: res(3) = [real(dp) :: -2, -2, 3]
        real(dp) :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_dp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_dp)
        print *, 'state = ',state%print()


    end subroutine test_dsolve
    subroutine test_qsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(qp) :: A(3,3) = transpose(reshape([real(qp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(qp) :: b  (3) = [real(qp) :: 1, 4, -1]
        real(qp) :: res(3) = [real(qp) :: -2, -2, 3]
        real(qp) :: x(3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_qp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_qp)
        print *, 'state = ',state%print()


    end subroutine test_qsolve

    !> Simple linear system with multiple right hand sides
    subroutine test_ssolve_multiple(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(sp) :: A(3,3) = transpose(reshape([real(sp) :: 1,-1, 2, &
                                                            0, 1, 1, &
                                                            1,-1, 3], [3,3]))
        real(sp) :: b(3,3) = transpose(reshape([real(sp) :: 0, 1, 2, &
                                                            1,-2,-1, &
                                                            2, 3,-1], [3,3]))
        real(sp) :: res(3,3) = transpose(reshape([real(sp) ::-5,-7,10, &
                                                           -1,-4, 2, &
                                                            2, 2,-3], [3,3]))
        real(sp) :: x(3,3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_sp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_sp)
        print *, 'state = ',state%print()


    end subroutine test_ssolve_multiple
    subroutine test_dsolve_multiple(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(dp) :: A(3,3) = transpose(reshape([real(dp) :: 1,-1, 2, &
                                                            0, 1, 1, &
                                                            1,-1, 3], [3,3]))
        real(dp) :: b(3,3) = transpose(reshape([real(dp) :: 0, 1, 2, &
                                                            1,-2,-1, &
                                                            2, 3,-1], [3,3]))
        real(dp) :: res(3,3) = transpose(reshape([real(dp) ::-5,-7,10, &
                                                           -1,-4, 2, &
                                                            2, 2,-3], [3,3]))
        real(dp) :: x(3,3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_dp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_dp)
        print *, 'state = ',state%print()


    end subroutine test_dsolve_multiple
    subroutine test_qsolve_multiple(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        real(qp) :: A(3,3) = transpose(reshape([real(qp) :: 1,-1, 2, &
                                                            0, 1, 1, &
                                                            1,-1, 3], [3,3]))
        real(qp) :: b(3,3) = transpose(reshape([real(qp) :: 0, 1, 2, &
                                                            1,-2,-1, &
                                                            2, 3,-1], [3,3]))
        real(qp) :: res(3,3) = transpose(reshape([real(qp) ::-5,-7,10, &
                                                           -1,-4, 2, &
                                                            2, 2,-3], [3,3]))
        real(qp) :: x(3,3)

        x = solve(a,b,err=state)
        error = state%error() .or. .not.all(abs(x-res)<abs(res*epsilon(0.0_qp)))

        print *, 'res = ',res
        print *, 'x   = ',x
        print *, 'err = ',abs(x-res)
        print *, 'tst = ',res*epsilon(0.0_qp)
        print *, 'state = ',state%print()


    end subroutine test_qsolve_multiple

    !> Complex linear system
    !> Militaru, Popa, "On the numerical solving of complex linear systems",
    !> Int J Pure Appl Math 76(1), 113-122, 2012.
    subroutine test_csolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        complex(sp) :: A(5,5), b(5), res(5), x(5)
        integer(ilp) :: i

        ! Fill in linear system
        A = (0.0_sp,0.0_sp)

        A(1:2,1) = [(19.73_sp,0.0_sp),(0.0_sp,-0.51_sp)]
        A(1:3,2) = [(12.11_sp,-1.0_sp),(32.3_sp,7.0_sp),(0.0_sp,-0.51_sp)]
        A(1:4,3) = [(0.0_sp,5.0_sp),(23.07_sp,0.0_sp),(70.0_sp,7.3_sp),(1.0_sp,1.1_sp)]
        A(2:5,4) = [(0.0_sp,1.0_sp),(3.95_sp,0.0_sp),(50.17_sp,0.0_sp),(0.0_sp,-9.351_sp)]
        A(3:5,5) = [(19.0_sp,31.83_sp),(45.51_sp,0.0_sp),(55.0_sp,0.0_sp)]

        b = [(77.38_sp,8.82_sp),(157.48_sp,19.8_sp),(1175.62_sp,20.69_sp),(912.12_sp,-801.75_sp),(550.0_sp,-1060.4_sp)]

        ! Exact result
        res = [(3.3_sp,-1.0_sp),(1.0_sp,0.17_sp),(5.5_sp,0.0_sp),(9.0_sp,0.0_sp),(10.0_sp,-17.75_sp)]

        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(abs(x-res)<abs(res)*1.0e-3_sp)

        do i=1,5
           print *, 'res = ',res(i),' x  =',x(i),' b =',b(i)
        end do
        print *, 'state = ',state%print()

    end subroutine test_csolve
    subroutine test_zsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        complex(dp) :: A(5,5), b(5), res(5), x(5)
        integer(ilp) :: i

        ! Fill in linear system
        A = (0.0_dp,0.0_dp)

        A(1:2,1) = [(19.73_dp,0.0_dp),(0.0_dp,-0.51_dp)]
        A(1:3,2) = [(12.11_dp,-1.0_dp),(32.3_dp,7.0_dp),(0.0_dp,-0.51_dp)]
        A(1:4,3) = [(0.0_dp,5.0_dp),(23.07_dp,0.0_dp),(70.0_dp,7.3_dp),(1.0_dp,1.1_dp)]
        A(2:5,4) = [(0.0_dp,1.0_dp),(3.95_dp,0.0_dp),(50.17_dp,0.0_dp),(0.0_dp,-9.351_dp)]
        A(3:5,5) = [(19.0_dp,31.83_dp),(45.51_dp,0.0_dp),(55.0_dp,0.0_dp)]

        b = [(77.38_dp,8.82_dp),(157.48_dp,19.8_dp),(1175.62_dp,20.69_dp),(912.12_dp,-801.75_dp),(550.0_dp,-1060.4_dp)]

        ! Exact result
        res = [(3.3_dp,-1.0_dp),(1.0_dp,0.17_dp),(5.5_dp,0.0_dp),(9.0_dp,0.0_dp),(10.0_dp,-17.75_dp)]

        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(abs(x-res)<abs(res)*1.0e-3_dp)

        do i=1,5
           print *, 'res = ',res(i),' x  =',x(i),' b =',b(i)
        end do
        print *, 'state = ',state%print()

    end subroutine test_zsolve
    subroutine test_wsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        complex(qp) :: A(5,5), b(5), res(5), x(5)
        integer(ilp) :: i

        ! Fill in linear system
        A = (0.0_qp,0.0_qp)

        A(1:2,1) = [(19.73_qp,0.0_qp),(0.0_qp,-0.51_qp)]
        A(1:3,2) = [(12.11_qp,-1.0_qp),(32.3_qp,7.0_qp),(0.0_qp,-0.51_qp)]
        A(1:4,3) = [(0.0_qp,5.0_qp),(23.07_qp,0.0_qp),(70.0_qp,7.3_qp),(1.0_qp,1.1_qp)]
        A(2:5,4) = [(0.0_qp,1.0_qp),(3.95_qp,0.0_qp),(50.17_qp,0.0_qp),(0.0_qp,-9.351_qp)]
        A(3:5,5) = [(19.0_qp,31.83_qp),(45.51_qp,0.0_qp),(55.0_qp,0.0_qp)]

        b = [(77.38_qp,8.82_qp),(157.48_qp,19.8_qp),(1175.62_qp,20.69_qp),(912.12_qp,-801.75_qp),(550.0_qp,-1060.4_qp)]

        ! Exact result
        res = [(3.3_qp,-1.0_qp),(1.0_qp,0.17_qp),(5.5_qp,0.0_qp),(9.0_qp,0.0_qp),(10.0_qp,-17.75_qp)]

        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(abs(x-res)<abs(res)*1.0e-3_qp)

        do i=1,5
           print *, 'res = ',res(i),' x  =',x(i),' b =',b(i)
        end do
        print *, 'state = ',state%print()

    end subroutine test_wsolve

    !> 2x2 Complex linear system
    !> https://math.stackexchange.com/questions/1996540/solving-linear-equation-systems-with-complex-coefficients-and-variables
    subroutine test_2x2_csolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        complex(sp) :: A(2,2), b(2), res(2), x(2)
        integer(ilp) :: i

        ! Fill in linear system
        A(1,:) = [(+1.0_sp,+1.0_sp),(-1.0_sp,0.0_sp)]
        A(2,:) = [(+1.0_sp,-1.0_sp),(+1.0_sp,1.0_sp)]

        b = [(0.0_sp,1.0_sp),(1.0_sp,0.0_sp)]

        ! Exact result
        res = [(0.5_sp,0.5_sp),(0.0_sp,0.0_sp)]

        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(abs(x-res)<max(tiny(0.0_sp),abs(res)*epsilon(0.0_sp)))

        do i=1,2
           print *, 'res = ',res(i),' x  =',x(i),' b =',b(i)
        end do
        print *, 'state = ',state%print()

    end subroutine test_2x2_csolve
    subroutine test_2x2_zsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        complex(dp) :: A(2,2), b(2), res(2), x(2)
        integer(ilp) :: i

        ! Fill in linear system
        A(1,:) = [(+1.0_dp,+1.0_dp),(-1.0_dp,0.0_dp)]
        A(2,:) = [(+1.0_dp,-1.0_dp),(+1.0_dp,1.0_dp)]

        b = [(0.0_dp,1.0_dp),(1.0_dp,0.0_dp)]

        ! Exact result
        res = [(0.5_dp,0.5_dp),(0.0_dp,0.0_dp)]

        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(abs(x-res)<max(tiny(0.0_dp),abs(res)*epsilon(0.0_dp)))

        do i=1,2
           print *, 'res = ',res(i),' x  =',x(i),' b =',b(i)
        end do
        print *, 'state = ',state%print()

    end subroutine test_2x2_zsolve
    subroutine test_2x2_wsolve(error)
        logical, intent(out) :: error

        type(linalg_state) :: state

        complex(qp) :: A(2,2), b(2), res(2), x(2)
        integer(ilp) :: i

        ! Fill in linear system
        A(1,:) = [(+1.0_qp,+1.0_qp),(-1.0_qp,0.0_qp)]
        A(2,:) = [(+1.0_qp,-1.0_qp),(+1.0_qp,1.0_qp)]

        b = [(0.0_qp,1.0_qp),(1.0_qp,0.0_qp)]

        ! Exact result
        res = [(0.5_qp,0.5_qp),(0.0_qp,0.0_qp)]

        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(abs(x-res)<max(tiny(0.0_qp),abs(res)*epsilon(0.0_qp)))

        do i=1,2
           print *, 'res = ',res(i),' x  =',x(i),' b =',b(i)
        end do
        print *, 'state = ',state%print()

    end subroutine test_2x2_wsolve

end module test_linalg_solve


