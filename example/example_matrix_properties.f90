! Example demonstrating matrix_rank, cond, slogdet and matrix_power
program example_matrix_properties
  use linear_algebra, only: dp, ilp, matrix_rank, cond, slogdet, matrix_power
  implicit none(type,external)
  real(dp) :: A(3,3), B(3,3), H(200,200), sgn, logabsdet
  integer(ilp) :: i

  A = reshape([4, 1, 0, &
               1, 3, 1, &
               0, 1, 2], [3,3])

  ! Rank-deficient: the third row is the sum of the first two
  B = A
  B(3,:) = B(1,:) + B(2,:)
  print '("rank(A) = ",i0,", rank(B) = ",i0)', matrix_rank(A), matrix_rank(B)

  ! Condition numbers in the 2-, 1- and Frobenius norms
  print '("cond(A) = ",3(f8.4,1x))', cond(A), cond(A,1_ilp), cond(A,'fro')

  ! A**3 and A**(-1)
  print '("A**3   = ",3(f9.4,1x))', transpose(matrix_power(A,3_ilp))
  print '("A**(-1)= ",3(f9.4,1x))', transpose(matrix_power(A,-1_ilp))

  ! An upper bidiagonal matrix with 0.01 on its diagonal: det(H) = 1e-400 underflows,
  ! its logarithm 200*log(0.01) = -921.03 does not
  H = 0
  do i = 1, 200
     H(i,i) = 0.01_dp
     if (i < 200) H(i,i+1) = 1
  end do
  call slogdet(H,sgn,logabsdet)
  print '("sign = ",f4.1,", log|det(H)| = ",f10.4)', sgn, logabsdet

end program example_matrix_properties
