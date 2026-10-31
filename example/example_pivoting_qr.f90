! Example demonstrating the QR factorization with column pivoting
program example_pivoting_qr
  use linear_algebra, only: dp, ilp, qr
  implicit none(type,external)
  real(dp) :: A(6,4), Q(6,4), R(4,4)
  integer(ilp) :: pivots(4), j

  ! A random matrix whose third column is twice the first: its rank is 3
  call random_number(A)
  A(:,3) = 2*A(:,1)

  ! Compute its pivoted QR factorization (reduced): A(:,pivots) = Q*R
  call qr(A,Q,R,pivots)

  ! Test factorization
  print *, maxval(abs(matmul(Q,R)-A(:,pivots)))

  ! |diag(R)| does not increase, and its last entry reveals the rank deficiency
  print *, 'pivots    =', pivots
  print *, '|diag(R)| =', [(abs(R(j,j)), j=1,4)]

end program example_pivoting_qr
