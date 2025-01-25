program example_getrf
  use linear_algebra, only: eye,dp,ilp,getrf
  implicit none(type,external)
  real(dp) :: A(3, 3)
  integer(ilp) :: ipiv(3),info
  
  A = eye(3)
  
  ! LAPACK matrix factorization interface (overwrite result)
  call getrf(size(A,1),size(A,2),A,size(A,1),ipiv,info)
  print *, info ! info==0: Success!

end program example_getrf
