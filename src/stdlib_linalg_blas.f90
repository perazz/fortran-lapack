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



          ! DOT forms the dot product of two vectors.
          ! uses unrolled loops for increments equal to one.
          interface dot
#ifdef STDLIB_EXTERNAL_BLAS
               real(dp) function ddot(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    integer(ilp) :: incx,incy,n,i 
                    real(dp) :: dx(*),dy(*) 
               end function ddot
#else
               module procedure stdlib_ddot
#endif
#ifdef STDLIB_EXTERNAL_BLAS
               real(sp) function sdot(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk 
                    implicit none(type,external) 
                    integer(ilp) :: incx,incy,n,i 
                    real(sp) :: sx(*),sy(*) 
               end function sdot
#else
               module procedure stdlib_sdot
#endif
          end interface dot



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





end module stdlib_linalg_blas
