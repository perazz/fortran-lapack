module test_linalg_solve
    use stdlib_linalg_constants
    use stdlib_linalg_state
    use stdlib_linalg_solve

    implicit none (type,external)




    contains

    !> Simple linear system
    subroutine test_solve(error)

        logical, intent(out) :: error

        type(linalg_state) :: state

        real(sp) :: A(3,3) = transpose(reshape([real(sp) :: 1, 3, 3, &
                                                            1, 3, 4, &
                                                            1, 4, 3], [3,3]))
        real(sp) :: b  (3) = [real(sp) :: 1, 4, -1]
        real(sp) :: res(3) = [real(sp) :: -2, -2, 3]
        real(sp) :: x(3)


        x = solve(a,b,err=state)

        error = state%error() .or. .not.all(x==res)

       print *, x






!>>> import imsl.linalg as la
!>>> a = [[1.0, 3.0, 3.0], [1.0, 3.0, 4.0], [1.0, 4.0, 3.0]]
!>>> b = [1.0, 4.0, -1.0]
!>>> x = la.lu_solve(a, b)
!>>> print(x)
![-2. -2. 3.]


    end subroutine test_solve



end module test_linalg_solve


