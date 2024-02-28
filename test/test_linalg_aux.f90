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

        type(linalg_state) :: state

        integer :: i
        integer, allocatable :: ijk(:,:),s(:)
        real, allocatable :: aa(:,:,:),x
        real :: a(5,10),b(6,10)
        type(array_descriptor) :: cfi1,cfi2

        allocate(aa(5,10,20))

        !s = strides(a(1:5:3,1:10:6,1:20:9))

        cfi1 = array_descriptor(a(1:5:3,1:10:6))
        call CFI_print(cfi1)
        cfi2 = array_descriptor(b(1:6:3,1:10:5))
        call CFI_print(cfi2)
stop

        error = .false.

        error = .not. size(strides(x))==0
        if (error) return

        error = .not. all(strides(a)==[1,1,1])
        if (error) return

        error = .not. all(strides(a)==strides(aa))
        if (error) return

!        do i=1,4
!            s = strides(a(::i,::2*i,::3*i))
!            error = .not. all(s==[1,2,3]*i)
!            print *, 'i=',i,' s=',s
!            if (error) return
!        end do









    end subroutine test_array_strides



end module test_linalg_aux


