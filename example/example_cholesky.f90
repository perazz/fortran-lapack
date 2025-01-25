! Example demonstrating Cholesky factorization as a subroutine
program example_cholesky
  use linear_algebra, only: cholesky
  implicit none

  real, dimension(3,3) :: A,L,U 

  ! Set real matrix
  A = reshape( [ [6, 15, 55], &
                 [15, 55, 225], &
                 [55, 225, 979] ], [3,3] )

  ! Decompose (lower)
  call cholesky(A, L, lower=.true.)

  ! Compare decomposition 
  print *, maxval(abs(A-matmul(L,transpose(L))))

  ! Decompose (upper)
  call cholesky(A, U, lower=.false.)
  
  ! Compare decomposition 
  print *, maxval(abs(A-matmul(transpose(U),U)))

end program example_cholesky

