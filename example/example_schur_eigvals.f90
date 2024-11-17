! This example includes eigenvalue computation in addition to 
! the Schur decomposition for a randomly generated matrix.
program example_schur_eigenvalues
  use stdlib_linalg_interface, only: schur, dp
  implicit none
  real(dp), allocatable :: A(:,:), T(:,:), Z(:,:)
  complex(dp), allocatable :: eigenvalues(:)
  integer :: n

  ! Create a random real-valued square matrix
  n = 5
  allocate(A(n,n), T(n,n), Z(n,n), eigenvalues(n))
  call random_number(A)

  ! Compute the Schur decomposition and eigenvalues
  call schur(A, T, Z, eigenvalues)

  ! Output results
  print *, "Random Matrix A:"
  print *, A
  print *, "Schur Form Matrix T:"
  print *, T
  print *, "Orthogonal Matrix Z:"
  print *, Z
  print *, "Eigenvalues:"
  print *, eigenvalues

  ! Test factorization: Z*T*Z^T = A
  print *, "Max error in reconstruction:", maxval(abs(matmul(Z, matmul(T, transpose(Z))) - A))

  deallocate(A, T, Z, eigenvalues)
end program example_schur_eigenvalues

