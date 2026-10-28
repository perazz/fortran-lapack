module test_la_aux
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

#ifdef LA_WITH_QP
        state = la_state(LINALG_SUCCESS,'128-bit real: ',1.0_qp)
        if (state%message/='128-bit real: 1.00000000000000000000000000000000000E+0000') error = .true.
#endif

        state = la_state(LINALG_SUCCESS,' 32-bit complex: ',(1.0_sp,1.0_sp))
        if (state%message/=' 32-bit complex: (1.00000000E+00,1.00000000E+00)') error = .true.

        state = la_state(LINALG_SUCCESS,' 64-bit complex: ',(1.0_dp,1.0_dp))
        if (state%message/=' 64-bit complex: (1.0000000000000000E+000,1.0000000000000000E+000)') error = .true.

#ifdef LA_WITH_QP
        state = la_state(LINALG_SUCCESS,'128-bit complex: ',(1.0_qp,1.0_qp))
        if (state%message/= &
        '128-bit complex: (1.00000000000000000000000000000000000E+0000,1.00000000000000000000000000000000000E+0000)') &
        error = .true.
#endif

        state = la_state(LINALG_SUCCESS,' 32-bit array: ',[(1.0_sp,0.0_sp),(0.0_sp,1.0_sp)])
        if (state%message/=' 32-bit array: [(1.00000000E+00,0.00000000E+00) (0.00000000E+00,1.00000000E+00)]') &
        error = .true.

        !> Arrays of any rank: their elements in array element order, then their shape
        state = la_state(LINALG_SUCCESS,'rank 2:',reshape([1,2,3,4,5,6],[2,3]))
        if (state%message/='rank 2: [1 2 3 4 5 6] (2x3)') error = .true.

        state = la_state(LINALG_SUCCESS,'rank 3:',reshape([1,2,3,4,5,6,7,8],[2,2,2]))
        if (state%message/='rank 3: [1 2 3 4 5 6 7 8] (2x2x2)') error = .true.

        state = la_state(LINALG_SUCCESS,'rank 14:',reshape([1,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,2]))
        if (state%message/='rank 14: [1 2] (1x1x1x1x1x1x1x1x1x1x1x1x1x2)') error = .true.

        state = la_state(LINALG_SUCCESS,'complex matrix:',reshape([(1.0_sp,0.0_sp),(0.0_sp,1.0_sp)],[1,2]))
        if (state%message/='complex matrix: [(1.00000000E+00,0.00000000E+00) (0.00000000E+00,1.00000000E+00)] (1x2)') &
        error = .true.

        !> Scalars and arrays of different ranks, in any order
        state = la_state(LINALG_SUCCESS,'n=',3,'a=',reshape([1,2,3,4],[2,2]),'v=',[5,6],'done')
        if (state%message/='n= 3 a= [1 2 3 4] (2x2) v= [5 6] done') error = .true.

        !> State flag with location
        state = la_state('test_formats',LINALG_SUCCESS,' 32-bit real: ',1.0_sp)
        if (state%print()/='[test_formats] returned Success!') error = .true.

    end subroutine test_formats

end module test_la_aux


