module test_linalg_aux
    use stdlib_linalg_interface
    use stdlib_linalg_arrays
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

        state = linalg_state(LINALG_SUCCESS,' 32-bit array: ',v1=[(1.0_sp,0.0_sp),(0.0_sp,1.0_sp)])
        if (state%message/=' 32-bit array: [(1.00000000E+00,0.00000000E+00) (0.00000000E+00,1.00000000E+00)]') &
        error = .true.

        !> State flag with location
        state = linalg_state('test_formats',LINALG_SUCCESS,' 32-bit real: ',1.0_sp)
        if (state%print()/='[test_formats] returned Success!') error = .true.

    end subroutine test_formats

    !> Test array strides
    subroutine test_array_strides(error)

        logical, intent(out) :: error

        integer :: i
        integer, allocatable :: ijk(:)
        integer :: m(10,10)

        error = .false.

        allocate(ijk(-21:349),source=0)

        error = .not. stride(ijk)==1
        if (error) return

        do i=1,340
            error = .not. stride(ijk(::i))==i
            if (error) return
            error = .not. stride(ijk(::-i))==-i
            if (error) return
        end do


        call rrr(m)
        call rrr(m(::4,::3))

    end subroutine test_array_strides

    subroutine rrr(mat)
        integer, intent(inout), target :: mat(:,:)
        integer, pointer :: p(:,:)

        type(array_descriptor) :: cfi_m,cfi_p

        integer :: n,m

        n = size(mat,1)
        m = size(mat,2)

        print *, 'shape(mat)=',shape(mat),' leading dim=',leading_dim(mat)


        p => mat

        cfi_m = array_descriptor(mat)
        cfi_p = array_descriptor(p)

        call CFI_print(cfi_m)
        call CFI_print(cfi_p)

    end subroutine rrr



end module test_linalg_aux


