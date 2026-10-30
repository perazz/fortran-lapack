! Example demonstrating a generalized least-squares fit with correlated errors
program example_generalized_lstsq
  use linear_algebra, only: generalized_lstsq, solve_generalized_lstsq
  implicit none

  integer, parameter :: m = 3
  real, dimension(m,2) :: A
  real, dimension(m,m) :: W
  real :: b(m), x(2)

  ! Straight line through 3 observations: intercept and slope
  A(:,1) = 1.0
  A(:,2) = [1.0, 2.0, 3.0]
  b = [1.0, 2.1, 2.9]

  ! Covariance matrix of the observation errors
  W = reshape([1.0,  0.5,  0.25, &
               0.5,  1.0,  0.5,  &
               0.25, 0.5,  1.0], [m,m])

  ! Function form
  x = generalized_lstsq(W, A, b)
  print '("intercept = ",f8.4,", slope = ",f8.4)', x

  ! Subroutine form, into an existing array
  call solve_generalized_lstsq(W, A, b, x)
  print '("intercept = ",f8.4,", slope = ",f8.4)', x

end program example_generalized_lstsq
