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
                    complex(sp) :: ca,cx(*),cy(*)
                    integer(ilp) :: incx,incy,n
               end subroutine caxpy
#else
               module procedure stdlib_caxpy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine daxpy(n,da,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: da,dx(*),dy(*)
                    integer(ilp) :: incx,incy,n
               end subroutine daxpy
#else
               module procedure stdlib_daxpy
#endif
               module procedure stdlib_qaxpy
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine saxpy(n,sa,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: sa,sx(*),sy(*)
                    integer(ilp) :: incx,incy,n
               end subroutine saxpy
#else
               module procedure stdlib_saxpy
#endif
               module procedure stdlib_waxpy
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zaxpy(n,za,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: za,zx(*),zy(*)
                    integer(ilp) :: incx,incy,n
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
                    integer(ilp) :: incx,incy,n
                    complex(sp) :: cx(*),cy(*)
               end subroutine ccopy
#else
               module procedure stdlib_ccopy
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dcopy(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(dp) :: dx(*),dy(*)
               end subroutine dcopy
#else
               module procedure stdlib_dcopy
#endif
               module procedure stdlib_qcopy
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine scopy(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(sp) :: sx(*),sy(*)
               end subroutine scopy
#else
               module procedure stdlib_scopy
#endif
               module procedure stdlib_wcopy
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zcopy(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(dp) :: zx(*),zy(*)
               end subroutine zcopy
#else
               module procedure stdlib_zcopy
#endif
          end interface copy

          ! DOT forms the dot product of two vectors.
          ! uses unrolled loops for increments equal to one.
          interface dot
#ifdef STDLIB_EXTERNAL_BLAS
               real(dp) function ddot(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(dp) :: dx(*),dy(*)
               end function ddot
#else
               module procedure stdlib_ddot
#endif
               module procedure stdlib_qdot
#ifdef STDLIB_EXTERNAL_BLAS
               real(sp) function sdot(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(sp) :: sx(*),sy(*)
               end function sdot
#else
               module procedure stdlib_sdot
#endif
          end interface dot

          ! DOTC forms the dot product of two complex vectors
          ! DOTC = X^H * Y
          interface dotc
#ifdef STDLIB_EXTERNAL_BLAS
               complex(sp) function cdotc(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(sp) :: cx(*),cy(*)
               end function cdotc
#else
               module procedure stdlib_cdotc
#endif
               module procedure stdlib_wdotc
#ifdef STDLIB_EXTERNAL_BLAS
               complex(dp) function zdotc(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(dp) :: zx(*),zy(*)
               end function zdotc
#else
               module procedure stdlib_zdotc
#endif
          end interface dotc

          ! DOTU forms the dot product of two complex vectors
          ! DOTU = X^T * Y
          interface dotu
#ifdef STDLIB_EXTERNAL_BLAS
               complex(sp) function cdotu(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(sp) :: cx(*),cy(*)
               end function cdotu
#else
               module procedure stdlib_cdotu
#endif
               module procedure stdlib_wdotu
#ifdef STDLIB_EXTERNAL_BLAS
               complex(dp) function zdotu(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(dp) :: zx(*),zy(*)
               end function zdotu
#else
               module procedure stdlib_zdotu
#endif
          end interface dotu

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
                    complex(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,kl,ku,lda,m,n
                    character :: trans
               end subroutine cgbmv
#else
               module procedure stdlib_cgbmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,kl,ku,lda,m,n
                    character :: trans
               end subroutine dgbmv
#else
               module procedure stdlib_dgbmv
#endif
               module procedure stdlib_qgbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,kl,ku,lda,m,n
                    character :: trans
               end subroutine sgbmv
#else
               module procedure stdlib_sgbmv
#endif
               module procedure stdlib_wgbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,kl,ku,lda,m,n
                    character :: trans
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
                    complex(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,m,n
                    character :: transa,transb
               end subroutine cgemm
#else
               module procedure stdlib_cgemm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,m,n
                    character :: transa,transb
               end subroutine dgemm
#else
               module procedure stdlib_dgemm
#endif
               module procedure stdlib_qgemm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,m,n
                    character :: transa,transb
               end subroutine sgemm
#else
               module procedure stdlib_sgemm
#endif
               module procedure stdlib_wgemm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,m,n
                    character :: transa,transb
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
                    complex(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
                    character :: trans
               end subroutine cgemv
#else
               module procedure stdlib_cgemv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
                    character :: trans
               end subroutine dgemv
#else
               module procedure stdlib_dgemv
#endif
               module procedure stdlib_qgemv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
                    character :: trans
               end subroutine sgemv
#else
               module procedure stdlib_sgemv
#endif
               module procedure stdlib_wgemv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
                    character :: trans
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
                    real(dp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
               end subroutine dger
#else
               module procedure stdlib_dger
#endif
               module procedure stdlib_qger
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sger(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
               end subroutine sger
#else
               module procedure stdlib_sger
#endif
          end interface ger

          ! GERC  performs the rank 1 operation
          ! A := alpha*x*y**H + A,
          ! where alpha is a scalar, x is an m element vector, y is an n element
          ! vector and A is an m by n matrix.
          interface gerc
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cgerc(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
               end subroutine cgerc
#else
               module procedure stdlib_cgerc
#endif
               module procedure stdlib_wgerc
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgerc(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
               end subroutine zgerc
#else
               module procedure stdlib_zgerc
#endif
          end interface gerc

          ! GERU  performs the rank 1 operation
          ! A := alpha*x*y**T + A,
          ! where alpha is a scalar, x is an m element vector, y is an n element
          ! vector and A is an m by n matrix.
          interface geru
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cgeru(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
               end subroutine cgeru
#else
               module procedure stdlib_cgeru
#endif
               module procedure stdlib_wgeru
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zgeru(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,m,n
               end subroutine zgeru
#else
               module procedure stdlib_zgeru
#endif
          end interface geru

          ! HBMV  performs the matrix-vector  operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n hermitian band matrix, with k super-diagonals.
          interface hbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine chbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,k,lda,n
                    character :: uplo
               end subroutine chbmv
#else
               module procedure stdlib_chbmv
#endif
               module procedure stdlib_whbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zhbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,k,lda,n
                    character :: uplo
               end subroutine zhbmv
#else
               module procedure stdlib_zhbmv
#endif
          end interface hbmv

          ! HEMM  performs one of the matrix-matrix operations
          ! C := alpha*A*B + beta*C,
          ! or
          ! C := alpha*B*A + beta*C,
          ! where alpha and beta are scalars, A is an hermitian matrix and  B and
          ! C are m by n matrices.
          interface hemm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine chemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: lda,ldb,ldc,m,n
                    character :: side,uplo
               end subroutine chemm
#else
               module procedure stdlib_chemm
#endif
               module procedure stdlib_whemm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: lda,ldb,ldc,m,n
                    character :: side,uplo
               end subroutine zhemm
#else
               module procedure stdlib_zhemm
#endif
          end interface hemm

          ! HEMV  performs the matrix-vector  operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n hermitian matrix.
          interface hemv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine chemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine chemv
#else
               module procedure stdlib_chemv
#endif
               module procedure stdlib_whemv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine zhemv
#else
               module procedure stdlib_zhemv
#endif
          end interface hemv

          ! HER   performs the hermitian rank 1 operation
          ! A := alpha*x*x**H + A,
          ! where alpha is a real scalar, x is an n element vector and A is an
          ! n by n hermitian matrix.
          interface her
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cher(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha
                    integer(ilp) :: incx,lda,n
                    character :: uplo
                    complex(sp) :: a(lda,*),x(*)
               end subroutine cher
#else
               module procedure stdlib_cher
#endif
               module procedure stdlib_wher
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zher(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha
                    integer(ilp) :: incx,lda,n
                    character :: uplo
                    complex(dp) :: a(lda,*),x(*)
               end subroutine zher
#else
               module procedure stdlib_zher
#endif
          end interface her

          ! HER2  performs the hermitian rank 2 operation
          ! A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
          ! where alpha is a scalar, x and y are n element vectors and A is an n
          ! by n hermitian matrix.
          interface her2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cher2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine cher2
#else
               module procedure stdlib_cher2
#endif
               module procedure stdlib_wher2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine zher2
#else
               module procedure stdlib_zher2
#endif
          end interface her2

          ! HER2K  performs one of the hermitian rank 2k operations
          ! C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
          ! or
          ! C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
          ! where  alpha and beta  are scalars with  beta  real,  C is an  n by n
          ! hermitian matrix and  A and B  are  n by k matrices in the first case
          ! and  k by n  matrices in the second case.
          interface her2k
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,a(lda,*),b(ldb,*),c(ldc,*)
                    real(sp) :: beta
                    integer(ilp) :: k,lda,ldb,ldc,n
                    character :: trans,uplo
               end subroutine cher2k
#else
               module procedure stdlib_cher2k
#endif
               module procedure stdlib_wher2k
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,a(lda,*),b(ldb,*),c(ldc,*)
                    real(dp) :: beta
                    integer(ilp) :: k,lda,ldb,ldc,n
                    character :: trans,uplo
               end subroutine zher2k
#else
               module procedure stdlib_zher2k
#endif
          end interface her2k

          ! HERK  performs one of the hermitian rank k operations
          ! C := alpha*A*A**H + beta*C,
          ! or
          ! C := alpha*A**H*A + beta*C,
          ! where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
          ! matrix and  A  is an  n by k  matrix in the  first case and a  k by n
          ! matrix in the second case.
          interface herk
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta
                    integer(ilp) :: k,lda,ldc,n
                    character :: trans,uplo
                    complex(sp) :: a(lda,*),c(ldc,*)
               end subroutine cherk
#else
               module procedure stdlib_cherk
#endif
               module procedure stdlib_wherk
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta
                    integer(ilp) :: k,lda,ldc,n
                    character :: trans,uplo
                    complex(dp) :: a(lda,*),c(ldc,*)
               end subroutine zherk
#else
               module procedure stdlib_zherk
#endif
          end interface herk

          ! HPMV  performs the matrix-vector operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n hermitian matrix, supplied in packed form.
          interface hpmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine chpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine chpmv
#else
               module procedure stdlib_chpmv
#endif
               module procedure stdlib_whpmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine zhpmv
#else
               module procedure stdlib_zhpmv
#endif
          end interface hpmv

          ! HPR    performs the hermitian rank 1 operation
          ! A := alpha*x*x**H + A,
          ! where alpha is a real scalar, x is an n element vector and A is an
          ! n by n hermitian matrix, supplied in packed form.
          interface hpr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine chpr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha
                    integer(ilp) :: incx,n
                    character :: uplo
                    complex(sp) :: ap(*),x(*)
               end subroutine chpr
#else
               module procedure stdlib_chpr
#endif
               module procedure stdlib_whpr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zhpr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha
                    integer(ilp) :: incx,n
                    character :: uplo
                    complex(dp) :: ap(*),x(*)
               end subroutine zhpr
#else
               module procedure stdlib_zhpr
#endif
          end interface hpr

          ! HPR2  performs the hermitian rank 2 operation
          ! A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
          ! where alpha is a scalar, x and y are n element vectors and A is an
          ! n by n hermitian matrix, supplied in packed form.
          interface hpr2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine chpr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine chpr2
#else
               module procedure stdlib_chpr2
#endif
               module procedure stdlib_whpr2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine zhpr2
#else
               module procedure stdlib_zhpr2
#endif
          end interface hpr2

          ! !
          ! NRM2 returns the euclidean norm of a vector via the function
          ! name, so that
          ! NRM2 := sqrt( x'*x )
          interface nrm2
#ifdef STDLIB_EXTERNAL_BLAS
               function dnrm2(n,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: dnrm2,x(*)
                    integer(ilp) :: incx,n
               end function dnrm2
#else
               module procedure stdlib_dnrm2
#endif
               module procedure stdlib_qnrm2
#ifdef STDLIB_EXTERNAL_BLAS
               function snrm2(n,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: snrm2,x(*)
                    integer(ilp) :: incx,n
               end function snrm2
#else
               module procedure stdlib_snrm2
#endif
          end interface nrm2

          ! ROT applies a plane rotation.
          interface rot
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine drot(n,dx,incx,dy,incy,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: c,s,dx(*),dy(*)
                    integer(ilp) :: incx,incy,n
               end subroutine drot
#else
               module procedure stdlib_drot
#endif
               module procedure stdlib_qrot
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine srot(n,sx,incx,sy,incy,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: c,s,sx(*),sy(*)
                    integer(ilp) :: incx,incy,n
               end subroutine srot
#else
               module procedure stdlib_srot
#endif
          end interface rot

          ! !
          ! The computation uses the formulas
          ! |x| = sqrt( Re(x)**2 + Im(x)**2 )
          ! sgn(x) = x / |x|  if x /= 0
          ! = 1        if x  = 0
          ! c = |a| / sqrt(|a|**2 + |b|**2)
          ! s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
          ! When a and b are real and r /= 0, the formulas simplify to
          ! r = sgn(a)*sqrt(|a|**2 + |b|**2)
          ! c = a / r
          ! s = b / r
          ! the same as in SROTG when |a| > |b|.  When |b| >= |a|, the
          ! sign of c and s will be different from those computed by SROTG
          ! if the signs of a and b are not the same.
          interface rotg
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine crotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: c
                    complex(sp) :: a,b,s
               end subroutine crotg
#else
               module procedure stdlib_crotg
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine drotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: a,b,c,s
               end subroutine drotg
#else
               module procedure stdlib_drotg
#endif
               module procedure stdlib_qrotg
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine srotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: a,b,c,s
               end subroutine srotg
#else
               module procedure stdlib_srotg
#endif
               module procedure stdlib_wrotg
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zrotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: c
                    complex(dp) :: a,b,s
               end subroutine zrotg
#else
               module procedure stdlib_zrotg
#endif
          end interface rotg

          ! APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
          ! (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
          ! (DY**T)
          ! DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
          ! LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
          ! WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
          ! DFLAG=-1._dp     DFLAG=0._dp        DFLAG=1._dp     DFLAG=-2.D0
          ! (DH11  DH12)    (1._dp  DH12)    (DH11  1._dp)    (1._dp  0._dp)
          ! H=(          )    (          )    (          )    (          )
          ! (DH21  DH22),   (DH21  1._dp),   (-1._dp DH22),   (0._dp  1._dp).
          ! SEE ROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
          interface rotm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine drotm(n,dx,incx,dy,incy,dparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(dp) :: dparam(5),dx(*),dy(*)
               end subroutine drotm
#else
               module procedure stdlib_drotm
#endif
               module procedure stdlib_qrotm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine srotm(n,sx,incx,sy,incy,sparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(sp) :: sparam(5),sx(*),sy(*)
               end subroutine srotm
#else
               module procedure stdlib_srotm
#endif
          end interface rotm

          ! CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
          ! THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(DD1)*DX1,SQRT(DD2)    DY2)**T.
          ! WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
          ! DFLAG=-1._dp     DFLAG=0._dp        DFLAG=1._dp     DFLAG=-2.D0
          ! (DH11  DH12)    (1._dp  DH12)    (DH11  1._dp)    (1._dp  0._dp)
          ! H=(          )    (          )    (          )    (          )
          ! (DH21  DH22),   (DH21  1._dp),   (-1._dp DH22),   (0._dp  1._dp).
          ! LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
          ! RESPECTIVELY. (VALUES OF 1._dp, -1._dp, OR 0._dp IMPLIED BY THE
          ! VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
          ! THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
          ! INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
          ! OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
          interface rotmg
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine drotmg(dd1,dd2,dx1,dy1,dparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: dd1,dd2,dx1,dy1,dparam(5)
               end subroutine drotmg
#else
               module procedure stdlib_drotmg
#endif
               module procedure stdlib_qrotmg
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine srotmg(sd1,sd2,sx1,sy1,sparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: sd1,sd2,sx1,sy1,sparam(5)
               end subroutine srotmg
#else
               module procedure stdlib_srotmg
#endif
          end interface rotmg

          ! SBMV  performs the matrix-vector  operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n symmetric band matrix, with k super-diagonals.
          interface sbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,k,lda,n
                    character :: uplo
               end subroutine dsbmv
#else
               module procedure stdlib_dsbmv
#endif
               module procedure stdlib_qsbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,k,lda,n
                    character :: uplo
               end subroutine ssbmv
#else
               module procedure stdlib_ssbmv
#endif
          end interface sbmv

          ! SCAL scales a vector by a constant.
          interface scal
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cscal(n,ca,cx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: ca,cx(*)
                    integer(ilp) :: incx,n
               end subroutine cscal
#else
               module procedure stdlib_cscal
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dscal(n,da,dx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: da,dx(*)
                    integer(ilp) :: incx,n
               end subroutine dscal
#else
               module procedure stdlib_dscal
#endif
               module procedure stdlib_qscal
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sscal(n,sa,sx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: sa,sx(*)
                    integer(ilp) :: incx,n
               end subroutine sscal
#else
               module procedure stdlib_sscal
#endif
               module procedure stdlib_wscal
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zscal(n,za,zx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: za,zx(*)
                    integer(ilp) :: incx,n
               end subroutine zscal
#else
               module procedure stdlib_zscal
#endif
          end interface scal

          ! Compute the inner product of two vectors with extended
          ! precision accumulation and result.
          ! Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
          ! SDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
          ! where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
          ! defined in a similar way using INCY.
          interface sdot
#ifdef STDLIB_EXTERNAL_BLAS
               real(dp) function dsdot(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(sp) :: sx(*),sy(*)
               end function dsdot
#else
               module procedure stdlib_dsdot
#endif
               module procedure stdlib_qsdot
          end interface sdot

          ! SPMV  performs the matrix-vector operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n symmetric matrix, supplied in packed form.
          interface spmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine dspmv
#else
               module procedure stdlib_dspmv
#endif
               module procedure stdlib_qspmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine sspmv
#else
               module procedure stdlib_sspmv
#endif
          end interface spmv

          ! SPR    performs the symmetric rank 1 operation
          ! A := alpha*x*x**T + A,
          ! where alpha is a real scalar, x is an n element vector and A is an
          ! n by n symmetric matrix, supplied in packed form.
          interface spr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dspr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,ap(*),x(*)
                    integer(ilp) :: incx,n
                    character :: uplo
               end subroutine dspr
#else
               module procedure stdlib_dspr
#endif
               module procedure stdlib_qspr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sspr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,ap(*),x(*)
                    integer(ilp) :: incx,n
                    character :: uplo
               end subroutine sspr
#else
               module procedure stdlib_sspr
#endif
          end interface spr

          ! SPR2  performs the symmetric rank 2 operation
          ! A := alpha*x*y**T + alpha*y*x**T + A,
          ! where alpha is a scalar, x and y are n element vectors and A is an
          ! n by n symmetric matrix, supplied in packed form.
          interface spr2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dspr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine dspr2
#else
               module procedure stdlib_dspr2
#endif
               module procedure stdlib_qspr2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sspr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,ap(*),x(*),y(*)
                    integer(ilp) :: incx,incy,n
                    character :: uplo
               end subroutine sspr2
#else
               module procedure stdlib_sspr2
#endif
          end interface spr2

          ! SROT applies a plane rotation, where the cos and sin (c and s) are real
          ! and the vectors cx and cy are complex.
          ! jack dongarra, linpack, 3/11/78.
          interface srot
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csrot(n,cx,incx,cy,incy,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(sp) :: c,s
                    complex(sp) :: cx,cy
               end subroutine csrot
#else
               module procedure stdlib_csrot
#endif
          end interface srot

          ! SSCAL scales a complex vector by a real constant.
          interface sscal
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csscal(n,sa,cx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: sa
                    integer(ilp) :: incx,n
                    complex(sp) :: cx(*)
               end subroutine csscal
#else
               module procedure stdlib_csscal
#endif
          end interface sscal

          ! SWAP interchanges two vectors.
          interface swap
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine cswap(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(sp) :: cx(*),cy(*)
               end subroutine cswap
#else
               module procedure stdlib_cswap
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dswap(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(dp) :: dx(*),dy(*)
               end subroutine dswap
#else
               module procedure stdlib_dswap
#endif
               module procedure stdlib_qswap
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine sswap(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    real(sp) :: sx(*),sy(*)
               end subroutine sswap
#else
               module procedure stdlib_sswap
#endif
               module procedure stdlib_wswap
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zswap(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,incy,n
                    complex(dp) :: zx(*),zy(*)
               end subroutine zswap
#else
               module procedure stdlib_zswap
#endif
          end interface swap

          ! SYMM  performs one of the matrix-matrix operations
          ! C := alpha*A*B + beta*C,
          ! or
          ! C := alpha*B*A + beta*C,
          ! where  alpha and beta are scalars, A is a symmetric matrix and  B and
          ! C are m by n matrices.
          interface symm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: lda,ldb,ldc,m,n
                    character :: side,uplo
               end subroutine csymm
#else
               module procedure stdlib_csymm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: lda,ldb,ldc,m,n
                    character :: side,uplo
               end subroutine dsymm
#else
               module procedure stdlib_dsymm
#endif
               module procedure stdlib_qsymm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: lda,ldb,ldc,m,n
                    character :: side,uplo
               end subroutine ssymm
#else
               module procedure stdlib_ssymm
#endif
               module procedure stdlib_wsymm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: lda,ldb,ldc,m,n
                    character :: side,uplo
               end subroutine zsymm
#else
               module procedure stdlib_zsymm
#endif
          end interface symm

          ! SYMV  performs the matrix-vector  operation
          ! y := alpha*A*x + beta*y,
          ! where alpha and beta are scalars, x and y are n element vectors and
          ! A is an n by n symmetric matrix.
          interface symv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine dsymv
#else
               module procedure stdlib_dsymv
#endif
               module procedure stdlib_qsymv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine ssymv
#else
               module procedure stdlib_ssymv
#endif
          end interface symv

          ! SYR   performs the symmetric rank 1 operation
          ! A := alpha*x*x**T + A,
          ! where alpha is a real scalar, x is an n element vector and A is an
          ! n by n symmetric matrix.
          interface syr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsyr(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,a(lda,*),x(*)
                    integer(ilp) :: incx,lda,n
                    character :: uplo
               end subroutine dsyr
#else
               module procedure stdlib_dsyr
#endif
               module procedure stdlib_qsyr
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssyr(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,a(lda,*),x(*)
                    integer(ilp) :: incx,lda,n
                    character :: uplo
               end subroutine ssyr
#else
               module procedure stdlib_ssyr
#endif
          end interface syr

          ! SYR2  performs the symmetric rank 2 operation
          ! A := alpha*x*y**T + alpha*y*x**T + A,
          ! where alpha is a scalar, x and y are n element vectors and A is an n
          ! by n symmetric matrix.
          interface syr2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine dsyr2
#else
               module procedure stdlib_dsyr2
#endif
               module procedure stdlib_qsyr2
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,a(lda,*),x(*),y(*)
                    integer(ilp) :: incx,incy,lda,n
                    character :: uplo
               end subroutine ssyr2
#else
               module procedure stdlib_ssyr2
#endif
          end interface syr2

          ! SYR2K  performs one of the symmetric rank 2k operations
          ! C := alpha*A*B**T + alpha*B*A**T + beta*C,
          ! or
          ! C := alpha*A**T*B + alpha*B**T*A + beta*C,
          ! where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
          ! and  A and B  are  n by k  matrices  in the  first  case  and  k by n
          ! matrices in the second case.
          interface syr2k
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,n
                    character :: trans,uplo
               end subroutine csyr2k
#else
               module procedure stdlib_csyr2k
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,n
                    character :: trans,uplo
               end subroutine dsyr2k
#else
               module procedure stdlib_dsyr2k
#endif
               module procedure stdlib_qsyr2k
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,n
                    character :: trans,uplo
               end subroutine ssyr2k
#else
               module procedure stdlib_ssyr2k
#endif
               module procedure stdlib_wsyr2k
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldb,ldc,n
                    character :: trans,uplo
               end subroutine zsyr2k
#else
               module procedure stdlib_zsyr2k
#endif
          end interface syr2k

          ! SYRK  performs one of the symmetric rank k operations
          ! C := alpha*A*A**T + beta*C,
          ! or
          ! C := alpha*A**T*A + beta*C,
          ! where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
          ! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
          ! in the second case.
          interface syrk
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine csyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,beta,a(lda,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldc,n
                    character :: trans,uplo
               end subroutine csyrk
#else
               module procedure stdlib_csyrk
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,beta,a(lda,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldc,n
                    character :: trans,uplo
               end subroutine dsyrk
#else
               module procedure stdlib_dsyrk
#endif
               module procedure stdlib_qsyrk
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ssyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,beta,a(lda,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldc,n
                    character :: trans,uplo
               end subroutine ssyrk
#else
               module procedure stdlib_ssyrk
#endif
               module procedure stdlib_wsyrk
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,beta,a(lda,*),c(ldc,*)
                    integer(ilp) :: k,lda,ldc,n
                    character :: trans,uplo
               end subroutine zsyrk
#else
               module procedure stdlib_zsyrk
#endif
          end interface syrk

          ! TBMV  performs one of the matrix-vector operations
          ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
          ! where x is an n element vector and  A is an n by n unit, or non-unit,
          ! upper or lower triangular band matrix, with ( k + 1 ) diagonals.
          interface tbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    complex(sp) :: a(lda,*),x(*)
               end subroutine ctbmv
#else
               module procedure stdlib_ctbmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    real(dp) :: a(lda,*),x(*)
               end subroutine dtbmv
#else
               module procedure stdlib_dtbmv
#endif
               module procedure stdlib_qtbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine stbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    real(sp) :: a(lda,*),x(*)
               end subroutine stbmv
#else
               module procedure stdlib_stbmv
#endif
               module procedure stdlib_wtbmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    complex(dp) :: a(lda,*),x(*)
               end subroutine ztbmv
#else
               module procedure stdlib_ztbmv
#endif
          end interface tbmv

          ! TBSV  solves one of the systems of equations
          ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
          ! where b and x are n element vectors and A is an n by n unit, or
          ! non-unit, upper or lower triangular band matrix, with ( k + 1 )
          ! diagonals.
          ! No test for singularity or near-singularity is included in this
          ! routine. Such tests must be performed before calling this routine.
          interface tbsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    complex(sp) :: a(lda,*),x(*)
               end subroutine ctbsv
#else
               module procedure stdlib_ctbsv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    real(dp) :: a(lda,*),x(*)
               end subroutine dtbsv
#else
               module procedure stdlib_dtbsv
#endif
               module procedure stdlib_qtbsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine stbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    real(sp) :: a(lda,*),x(*)
               end subroutine stbsv
#else
               module procedure stdlib_stbsv
#endif
               module procedure stdlib_wtbsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,k,lda,n
                    character :: diag,trans,uplo
                    complex(dp) :: a(lda,*),x(*)
               end subroutine ztbsv
#else
               module procedure stdlib_ztbsv
#endif
          end interface tbsv

          ! TPMV  performs one of the matrix-vector operations
          ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
          ! where x is an n element vector and  A is an n by n unit, or non-unit,
          ! upper or lower triangular matrix, supplied in packed form.
          interface tpmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    complex(sp) :: ap(*),x(*)
               end subroutine ctpmv
#else
               module procedure stdlib_ctpmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    real(dp) :: ap(*),x(*)
               end subroutine dtpmv
#else
               module procedure stdlib_dtpmv
#endif
               module procedure stdlib_qtpmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine stpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    real(sp) :: ap(*),x(*)
               end subroutine stpmv
#else
               module procedure stdlib_stpmv
#endif
               module procedure stdlib_wtpmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    complex(dp) :: ap(*),x(*)
               end subroutine ztpmv
#else
               module procedure stdlib_ztpmv
#endif
          end interface tpmv

          ! TPSV  solves one of the systems of equations
          ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
          ! where b and x are n element vectors and A is an n by n unit, or
          ! non-unit, upper or lower triangular matrix, supplied in packed form.
          ! No test for singularity or near-singularity is included in this
          ! routine. Such tests must be performed before calling this routine.
          interface tpsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    complex(sp) :: ap(*),x(*)
               end subroutine ctpsv
#else
               module procedure stdlib_ctpsv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    real(dp) :: ap(*),x(*)
               end subroutine dtpsv
#else
               module procedure stdlib_dtpsv
#endif
               module procedure stdlib_qtpsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine stpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    real(sp) :: ap(*),x(*)
               end subroutine stpsv
#else
               module procedure stdlib_stpsv
#endif
               module procedure stdlib_wtpsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,n
                    character :: diag,trans,uplo
                    complex(dp) :: ap(*),x(*)
               end subroutine ztpsv
#else
               module procedure stdlib_ztpsv
#endif
          end interface tpsv

          ! TRMM  performs one of the matrix-matrix operations
          ! B := alpha*op( A )*B,   or   B := alpha*B*op( A )
          ! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
          ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
          ! op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
          interface trmm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine ctrmm
#else
               module procedure stdlib_ctrmm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine dtrmm
#else
               module procedure stdlib_dtrmm
#endif
               module procedure stdlib_qtrmm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine strmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine strmm
#else
               module procedure stdlib_strmm
#endif
               module procedure stdlib_wtrmm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine ztrmm
#else
               module procedure stdlib_ztrmm
#endif
          end interface trmm

          ! TRMV  performs one of the matrix-vector operations
          ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
          ! where x is an n element vector and  A is an n by n unit, or non-unit,
          ! upper or lower triangular matrix.
          interface trmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctrmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    complex(sp) :: a(lda,*),x(*)
               end subroutine ctrmv
#else
               module procedure stdlib_ctrmv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtrmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    real(dp) :: a(lda,*),x(*)
               end subroutine dtrmv
#else
               module procedure stdlib_dtrmv
#endif
               module procedure stdlib_qtrmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine strmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    real(sp) :: a(lda,*),x(*)
               end subroutine strmv
#else
               module procedure stdlib_strmv
#endif
               module procedure stdlib_wtrmv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztrmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    complex(dp) :: a(lda,*),x(*)
               end subroutine ztrmv
#else
               module procedure stdlib_ztrmv
#endif
          end interface trmv

          ! TRSM  solves one of the matrix equations
          ! op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
          ! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
          ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
          ! op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
          ! The matrix X is overwritten on B.
          interface trsm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine ctrsm
#else
               module procedure stdlib_ctrsm
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine dtrsm
#else
               module procedure stdlib_dtrsm
#endif
               module procedure stdlib_qtrsm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine strsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine strsm
#else
               module procedure stdlib_strsm
#endif
               module procedure stdlib_wtrsm
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp) :: alpha,a(lda,*),b(ldb,*)
                    integer(ilp) :: lda,ldb,m,n
                    character :: diag,side,transa,uplo
               end subroutine ztrsm
#else
               module procedure stdlib_ztrsm
#endif
          end interface trsm

          ! TRSV  solves one of the systems of equations
          ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
          ! where b and x are n element vectors and A is an n by n unit, or
          ! non-unit, upper or lower triangular matrix.
          ! No test for singularity or near-singularity is included in this
          ! routine. Such tests must be performed before calling this routine.
          interface trsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ctrsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    complex(sp) :: a(lda,*),x(*)
               end subroutine ctrsv
#else
               module procedure stdlib_ctrsv
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    real(dp) :: a(lda,*),x(*)
               end subroutine dtrsv
#else
               module procedure stdlib_dtrsv
#endif
               module procedure stdlib_qtrsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine strsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    real(sp) :: a(lda,*),x(*)
               end subroutine strsv
#else
               module procedure stdlib_strsv
#endif
               module procedure stdlib_wtrsv
#ifdef STDLIB_EXTERNAL_BLAS
               subroutine ztrsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp) :: incx,lda,n
                    character :: diag,trans,uplo
                    complex(dp) :: a(lda,*),x(*)
               end subroutine ztrsv
#else
               module procedure stdlib_ztrsv
#endif
          end interface trsv

end module stdlib_linalg_blas
