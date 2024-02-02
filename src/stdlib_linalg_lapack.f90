module stdlib_linalg_lapack
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack_aux
     use stdlib_linalg_lapack_s
     use stdlib_linalg_lapack_d
     use stdlib_linalg_lapack_q
     use stdlib_linalg_lapack_c
     use stdlib_linalg_lapack_z
     use stdlib_linalg_lapack_w
     implicit none(type, external)
     public

          ! ROT   applies a plane rotation, where the cos (C) is real and the
          ! sin (S) is complex, and the vectors CX and CY are complex.
          interface rot
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine crot(n, cx, incx, cy, incy, c, s)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    integer(ilp) :: incx, incy, n
                    real(sp) :: c
                    complex(sp) :: s, cx, cy
               end subroutine crot
#else
               module procedure stdlib_crot
#endif
               module procedure stdlib_wrot
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zrot(n, cx, incx, cy, incy, c, s)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    integer(ilp) :: incx, incy, n
                    real(dp) :: c
                    complex(dp) :: s, cx, cy
               end subroutine zrot
#else
               module procedure stdlib_zrot
#endif
          end interface rot

          ! SPMV  performs the matrix-vector operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n symmetric matrix, supplied in packed form.
          interface spmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, incy, n
                    complex(sp) :: alpha, beta, ap, x, y
               end subroutine cspmv
#else
               module procedure stdlib_cspmv
#endif
               module procedure stdlib_wspmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, incy, n
                    complex(dp) :: alpha, beta, ap, x, y
               end subroutine zspmv
#else
               module procedure stdlib_zspmv
#endif
          end interface spmv

          ! SPR    performs the symmetric rank 1 operation
          ! A := alpha*x*x**H + A,
          ! where alpha is a complex scalar, x is an n element vector and A is an
          ! n by n symmetric matrix, supplied in packed form.
          interface spr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cspr(uplo, n, alpha, x, incx, ap)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, n
                    complex(sp) :: alpha, ap, x
               end subroutine cspr
#else
               module procedure stdlib_cspr
#endif
               module procedure stdlib_wspr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zspr(uplo, n, alpha, x, incx, ap)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, n
                    complex(dp) :: alpha, ap, x
               end subroutine zspr
#else
               module procedure stdlib_zspr
#endif
          end interface spr

          ! SYMV  performs the matrix-vector  operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n symmetric matrix.
          interface symv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, incy, lda, n
                    complex(sp) :: alpha, beta, a, x, y
               end subroutine csymv
#else
               module procedure stdlib_csymv
#endif
               module procedure stdlib_wsymv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, incy, lda, n
                    complex(dp) :: alpha, beta, a, x, y
               end subroutine zsymv
#else
               module procedure stdlib_zsymv
#endif
          end interface symv

          ! SYR   performs the symmetric rank 1 operation
          ! A := alpha*x*x**H + A,
          ! where alpha is a complex scalar, x is an n element vector and A is an
          ! n by n symmetric matrix.
          interface syr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csyr(uplo, n, alpha, x, incx, a, lda)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, lda, n
                    complex(sp) :: alpha, a, x
               end subroutine csyr
#else
               module procedure stdlib_csyr
#endif
               module procedure stdlib_wsyr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zsyr(uplo, n, alpha, x, incx, a, lda)
                    import sp, dp, qp, ilp, lk
                    implicit none(type, external)
                    character :: uplo
                    integer(ilp) :: incx, lda, n
                    complex(dp) :: alpha, a, x
               end subroutine zsyr
#else
               module procedure stdlib_zsyr
#endif
          end interface syr

end module stdlib_linalg_lapack
