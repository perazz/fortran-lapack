module test_linalg_aux
    use stdlib_linalg_state
    implicit none (type,external)




    contains


    !> Test output formats of real/complex numbers
    subroutine test_formats(error)

        logical, intent(out) :: error

        type(linalg_state) :: state

        error = .false.

        state = linalg_state(LINALG_SUCCESS,' 32-bit real: ',1.0_sp)
        if (state%message/=' 32-bit real: 1.00000000E+00') error = .true.

        state = linalg_state(LINALG_SUCCESS,' 64-bit real: ',1.0_dp)
        if (state%message/=' 64-bit real: 1.0000000000000000E+000') error = .true.

        state = linalg_state(LINALG_SUCCESS,'128-bit real: ',1.0_qp)
        if (state%message/='128-bit real: 1.00000000000000000000000000000000000E+0000') error = .true.

        state = linalg_state(LINALG_SUCCESS,' 32-bit complex: ',(1.0_sp,1.0_sp))
        if (state%message/=' 32-bit complex: (1.00000000E+00,1.00000000E+00)') error = .true.

        state = linalg_state(LINALG_SUCCESS,' 64-bit complex: ',(1.0_dp,1.0_dp))
        if (state%message/=' 64-bit complex: (1.0000000000000000E+000,1.0000000000000000E+000)') error = .true.

        state = linalg_state(LINALG_SUCCESS,'128-bit complex: ',(1.0_qp,1.0_qp))
        if (state%message/= &
        '128-bit complex: (1.00000000000000000000000000000000000E+0000,1.00000000000000000000000000000000000E+0000)') &
        error = .true.

    end subroutine test_formats



end module test_linalg_aux


