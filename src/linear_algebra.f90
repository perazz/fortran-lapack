module linear_algebra
     use la_blas
     use la_constants
     use la_cholesky, only: chol, cholesky
     use la_determinant
     use la_eye, only: eye, diag
     use la_inverse
     use la_lapack
     use la_least_squares
     use la_solve
     use la_state_type
     use la_svd
     use la_eig
     use la_qr
     use la_norms
     use la_schur
     use la_pseudoinverse
     implicit none(type,external)
     public
     
     !> Constructs the identity matrix.
     !! This interface provides procedures to generate an identity matrix of a given size.
     !! The resulting matrix has 1s on the diagonal and 0s elsewhere.     
     public :: eye

     
end module linear_algebra
