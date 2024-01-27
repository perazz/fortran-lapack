module stdlib_linalg_blas
     use stdlib_linalg_constants
     use stdlib_linalg_blas_aux
     use stdlib_linalg_blas_s
     use stdlib_linalg_blas_d
     use stdlib_linalg_blas_q
     use stdlib_linalg_blas_c
     use stdlib_linalg_blas_z
     use stdlib_linalg_blas_w
     implicit none(type,external)
     public

          ! GEMV performs one of the matrix-vector operations
          ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
          ! y := alpha*A**H*x + beta*y,
          ! where alpha and beta are scalars, x and y are vectors and A is an
          ! m by n matrix.
          interface gemv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(sp) :: alpha 
                    complex(sp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    complex(sp) :: a(lda,*) 
                    complex(sp) :: x(*) 
                    complex(sp) :: y(*) 
               end subroutine cgemv
#else
               module procedure stdlib_cgemv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(dp) :: alpha 
                    real(dp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    real(dp) :: a(lda,*) 
                    real(dp) :: x(*) 
                    real(dp) :: y(*) 
               end subroutine dgemv
#else
               module procedure stdlib_dgemv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(sp) :: alpha 
                    real(sp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    real(sp) :: a(lda,*) 
                    real(sp) :: x(*) 
                    real(sp) :: y(*) 
               end subroutine sgemv
#else
               module procedure stdlib_sgemv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(dp) :: alpha 
                    complex(dp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    complex(dp) :: a(lda,*) 
                    complex(dp) :: x(*) 
                    complex(dp) :: y(*) 
               end subroutine zgemv
#else
               module procedure stdlib_zgemv
#endif
          end interface gemv



end module stdlib_linalg_blas
