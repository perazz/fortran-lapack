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

          ! AXPY constant times a vector plus a vector.
          interface axpy
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine caxpy(n,ca,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(sp) :: ca 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    complex(sp) :: cx(*) 
                    complex(sp) :: cy(*) 
               end subroutine caxpy
#else
               module procedure stdlib_caxpy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine daxpy(n,da,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(dp) :: da 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    real(dp) :: dx(*) 
                    real(dp) :: dy(*) 
               end subroutine daxpy
#else
               module procedure stdlib_daxpy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine saxpy(n,sa,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(sp) :: sa 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    real(sp) :: sx(*) 
                    real(sp) :: sy(*) 
               end subroutine saxpy
#else
               module procedure stdlib_saxpy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zaxpy(n,za,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(dp) :: za 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    complex(dp) :: zx(*) 
                    complex(dp) :: zy(*) 
               end subroutine zaxpy
#else
               module procedure stdlib_zaxpy
#endif
          end interface axpy



          ! COPY copies a vector x to a vector y.
          interface copy
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ccopy(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    complex(sp) :: cx(*) 
                    complex(sp) :: cy(*) 
               end subroutine ccopy
#else
               module procedure stdlib_ccopy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dcopy(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    real(dp) :: dx(*) 
                    real(dp) :: dy(*) 
               end subroutine dcopy
#else
               module procedure stdlib_dcopy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine scopy(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    real(sp) :: sx(*) 
                    real(sp) :: sy(*) 
               end subroutine scopy
#else
               module procedure stdlib_scopy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zcopy(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: n 
                    complex(dp) :: zx(*) 
                    complex(dp) :: zy(*) 
               end subroutine zcopy
#else
               module procedure stdlib_zcopy
#endif
          end interface copy



          ! GBMV  performs one of the matrix-vector operations
          ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
          ! y := alpha*A**H*x + beta*y,
          ! where alpha and beta are scalars, x and y are vectors and A is an
          ! m by n band matrix, with kl sub-diagonals and ku super-diagonals.
          interface gbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(sp) :: alpha 
                    complex(sp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: kl 
                    integer(ilp) :: ku 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    complex(sp) :: a(lda,*) 
                    complex(sp) :: x(*) 
                    complex(sp) :: y(*) 
               end subroutine cgbmv
#else
               module procedure stdlib_cgbmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(dp) :: alpha 
                    real(dp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: kl 
                    integer(ilp) :: ku 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    real(dp) :: a(lda,*) 
                    real(dp) :: x(*) 
                    real(dp) :: y(*) 
               end subroutine dgbmv
#else
               module procedure stdlib_dgbmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(sp) :: alpha 
                    real(sp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: kl 
                    integer(ilp) :: ku 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    real(sp) :: a(lda,*) 
                    real(sp) :: x(*) 
                    real(sp) :: y(*) 
               end subroutine sgbmv
#else
               module procedure stdlib_sgbmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(dp) :: alpha 
                    complex(dp) :: beta 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: kl 
                    integer(ilp) :: ku 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: trans 
                    complex(dp) :: a(lda,*) 
                    complex(dp) :: x(*) 
                    complex(dp) :: y(*) 
               end subroutine zgbmv
#else
               module procedure stdlib_zgbmv
#endif
          end interface gbmv



          ! GEMM  performs one of the matrix-matrix operations
          ! C := alpha*op( A )*op( B ) + beta*C,
          ! where  op( X ) is one of
          ! op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
          ! alpha and beta are scalars, and A, B and C are matrices, with op( A )
          ! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
          interface gemm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(sp) :: alpha 
                    complex(sp) :: beta 
                    integer(ilp) :: k 
                    integer(ilp) :: lda 
                    integer(ilp) :: ldb 
                    integer(ilp) :: ldc 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: transa 
                    character :: transb 
                    complex(sp) :: a(lda,*) 
                    complex(sp) :: b(ldb,*) 
                    complex(sp) :: c(ldc,*) 
               end subroutine cgemm
#else
               module procedure stdlib_cgemm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(dp) :: alpha 
                    real(dp) :: beta 
                    integer(ilp) :: k 
                    integer(ilp) :: lda 
                    integer(ilp) :: ldb 
                    integer(ilp) :: ldc 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: transa 
                    character :: transb 
                    real(dp) :: a(lda,*) 
                    real(dp) :: b(ldb,*) 
                    real(dp) :: c(ldc,*) 
               end subroutine dgemm
#else
               module procedure stdlib_dgemm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(sp) :: alpha 
                    real(sp) :: beta 
                    integer(ilp) :: k 
                    integer(ilp) :: lda 
                    integer(ilp) :: ldb 
                    integer(ilp) :: ldc 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: transa 
                    character :: transb 
                    real(sp) :: a(lda,*) 
                    real(sp) :: b(ldb,*) 
                    real(sp) :: c(ldc,*) 
               end subroutine sgemm
#else
               module procedure stdlib_sgemm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    complex(dp) :: alpha 
                    complex(dp) :: beta 
                    integer(ilp) :: k 
                    integer(ilp) :: lda 
                    integer(ilp) :: ldb 
                    integer(ilp) :: ldc 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    character :: transa 
                    character :: transb 
                    complex(dp) :: a(lda,*) 
                    complex(dp) :: b(ldb,*) 
                    complex(dp) :: c(ldc,*) 
               end subroutine zgemm
#else
               module procedure stdlib_zgemm
#endif
          end interface gemm



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



          ! GER   performs the rank 1 operation
          ! A := alpha*x*y**T + A,
          ! where alpha is a scalar, x is an m element vector, y is an n element
          ! vector and A is an m by n matrix.
          interface ger
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(dp) :: alpha 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    real(dp) :: a(lda,*) 
                    real(dp) :: x(*) 
                    real(dp) :: y(*) 
               end subroutine dger
#else
               module procedure stdlib_dger
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sger(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    real(sp) :: alpha 
                    integer(ilp) :: incx 
                    integer(ilp) :: incy 
                    integer(ilp) :: lda 
                    integer(ilp) :: m 
                    integer(ilp) :: n 
                    real(sp) :: a(lda,*) 
                    real(sp) :: x(*) 
                    real(sp) :: y(*) 
               end subroutine sger
#else
               module procedure stdlib_sger
#endif
          end interface ger





end module stdlib_linalg_blas
