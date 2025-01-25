! QR example with pre-allocated storage
program example_qr2
  use la_interface, only: dp, ilp, qr, qr_space, linalg_state
  implicit none(type,external)
  real(dp) :: A(104, 32), Q(104,32), R(32,32)
  real(dp), allocatable :: work(:)
  integer(ilp) :: lwork
  type(linalg_state) :: err
  
  ! Create a random matrix
  call random_number(A)

  ! Prepare QR workspace
  call qr_space(A,lwork)
  allocate(work(lwork))

  ! Compute its QR factorization (reduced)
  call qr(A,Q,R,storage=work,err=err)

  ! Test factorization: Q*R = A 
  print *, maxval(abs(matmul(Q,R)-A))
  print *, err%print() 

end program example_qr2
