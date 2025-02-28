program example_qr
  use linear_algebra, only: dp, qr
  implicit none(type,external)
  real(dp) :: A(104, 32), Q(104,32), R(32,32)
  
  ! Create a random matrix
  call random_number(A)

  ! Compute its QR factorization (reduced)
  call qr(A,Q,R)

  ! Test factorization: Q*R = A 
  print *, maxval(abs(matmul(Q,R)-A)) 

end program example_qr
