module test_linalg_aux
    use linear_algebra
    implicit none (type,external)

    contains

    !> Test output formats of real/complex numbers
    subroutine test_formats(error)

        logical, intent(out) :: error

        type(la_state) :: state

        error = .false.

        state = la_state(LINALG_SUCCESS,' 32-bit real: ',1.0_sp)
        if (state%message/=' 32-bit real: 1.00000000E+00') error = .true.

        state = la_state(LINALG_SUCCESS,' 64-bit real: ',1.0_dp)
        if (state%message/=' 64-bit real: 1.0000000000000000E+000') error = .true.

        state = la_state(LINALG_SUCCESS,'128-bit real: ',1.0_qp)
        if (state%message/='128-bit real: 1.00000000000000000000000000000000000E+0000') error = .true.

        state = la_state(LINALG_SUCCESS,' 32-bit complex: ',(1.0_sp,1.0_sp))
        if (state%message/=' 32-bit complex: (1.00000000E+00,1.00000000E+00)') error = .true.

        state = la_state(LINALG_SUCCESS,' 64-bit complex: ',(1.0_dp,1.0_dp))
        if (state%message/=' 64-bit complex: (1.0000000000000000E+000,1.0000000000000000E+000)') error = .true.

        state = la_state(LINALG_SUCCESS,'128-bit complex: ',(1.0_qp,1.0_qp))
        if (state%message/= &
        '128-bit complex: (1.00000000000000000000000000000000000E+0000,1.00000000000000000000000000000000000E+0000)') &
        error = .true.

        state = la_state(LINALG_SUCCESS,' 32-bit array: ',[(1.0_sp,0.0_sp),(0.0_sp,1.0_sp)])
        if (state%message/=' 32-bit array: [(1.00000000E+00,0.00000000E+00) (0.00000000E+00,1.00000000E+00)]') &
        error = .true.

        !> State flag with location
        state = la_state('test_formats',LINALG_SUCCESS,' 32-bit real: ',1.0_sp)
        if (state%print()/='[test_formats] returned Success!') error = .true.

    end subroutine test_formats

end module test_linalg_aux


