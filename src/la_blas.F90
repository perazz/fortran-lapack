!> Precision-agnostic BLAS interface
module la_blas
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_ban
     use la_blas_level2_gen
     use la_blas_level2_pac
     use la_blas_level2_sym
     use la_blas_level2_tri
     use la_blas_level3_gen
     use la_blas_level3_sym
     use la_blas_level3_tri
     implicit none(type,external)
     public

          !> AXPY: constant times a vector plus a vector.
          interface axpy
#ifdef LA_EXTERNAL_BLAS
               pure subroutine caxpy(n,ca,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: ca,cx(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(sp),intent(inout) :: cy(*)
               end subroutine caxpy
#else
               module procedure la_caxpy
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine daxpy(n,da,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: da,dx(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(dp),intent(inout) :: dy(*)
               end subroutine daxpy
#else
               module procedure la_daxpy
#endif
#ifdef LA_WITH_QP
               module procedure la_qaxpy
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine saxpy(n,sa,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: sa,sx(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(inout) :: sy(*)
               end subroutine saxpy
#else
               module procedure la_saxpy
#endif
#ifdef LA_WITH_QP
               module procedure la_waxpy
#endif
#ifdef LA_WITH_XDP
               module procedure la_xaxpy
#endif
#ifdef LA_WITH_XDP
               module procedure la_yaxpy
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zaxpy(n,za,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: za,zx(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(dp),intent(inout) :: zy(*)
               end subroutine zaxpy
#else
               module procedure la_zaxpy
#endif
          end interface axpy

          !> COPY: copies a vector x to a vector y.
          interface copy
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ccopy(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(sp),intent(in) :: cx(*)
                    complex(sp),intent(out) :: cy(*)
               end subroutine ccopy
#else
               module procedure la_ccopy
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dcopy(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(dp),intent(in) :: dx(*)
                    real(dp),intent(out) :: dy(*)
               end subroutine dcopy
#else
               module procedure la_dcopy
#endif
#ifdef LA_WITH_QP
               module procedure la_qcopy
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine scopy(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(in) :: sx(*)
                    real(sp),intent(out) :: sy(*)
               end subroutine scopy
#else
               module procedure la_scopy
#endif
#ifdef LA_WITH_QP
               module procedure la_wcopy
#endif
#ifdef LA_WITH_XDP
               module procedure la_xcopy
#endif
#ifdef LA_WITH_XDP
               module procedure la_ycopy
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zcopy(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(dp),intent(in) :: zx(*)
                    complex(dp),intent(out) :: zy(*)
               end subroutine zcopy
#else
               module procedure la_zcopy
#endif
          end interface copy

          !> DOT: forms the dot product of two vectors.
          !> uses unrolled loops for increments equal to one.
          interface dot
#ifdef LA_EXTERNAL_BLAS
               pure real(dp) function ddot(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(dp),intent(in) :: dx(*),dy(*)
               end function ddot
#else
               module procedure la_ddot
#endif
#ifdef LA_WITH_QP
               module procedure la_qdot
#endif
#ifdef LA_EXTERNAL_BLAS
               pure real(sp) function sdot(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(in) :: sx(*),sy(*)
               end function sdot
#else
               module procedure la_sdot
#endif
#ifdef LA_WITH_XDP
               module procedure la_xdot
#endif
          end interface dot

          !> DOTC: forms the dot product of two complex vectors
          !> DOTC = X^H * Y
          interface dotc
#ifdef LA_EXTERNAL_BLAS
               pure complex(sp) function cdotc(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(sp),intent(in) :: cx(*),cy(*)
               end function cdotc
#else
               module procedure la_cdotc
#endif
#ifdef LA_WITH_QP
               module procedure la_wdotc
#endif
#ifdef LA_WITH_XDP
               module procedure la_ydotc
#endif
#ifdef LA_EXTERNAL_BLAS
               pure complex(dp) function zdotc(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(dp),intent(in) :: zx(*),zy(*)
               end function zdotc
#else
               module procedure la_zdotc
#endif
          end interface dotc

          !> DOTU: forms the dot product of two complex vectors
          !> DOTU = X^T * Y
          interface dotu
#ifdef LA_EXTERNAL_BLAS
               pure complex(sp) function cdotu(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(sp),intent(in) :: cx(*),cy(*)
               end function cdotu
#else
               module procedure la_cdotu
#endif
#ifdef LA_WITH_QP
               module procedure la_wdotu
#endif
#ifdef LA_WITH_XDP
               module procedure la_ydotu
#endif
#ifdef LA_EXTERNAL_BLAS
               pure complex(dp) function zdotu(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(dp),intent(in) :: zx(*),zy(*)
               end function zdotu
#else
               module procedure la_zdotu
#endif
          end interface dotu

          !> GBMV:  performs one of the matrix-vector operations
          !> y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
          !> y := alpha*A**H*x + beta*y,
          !> where alpha and beta are scalars, x and y are vectors and A is an
          !> m by n band matrix, with kl sub-diagonals and ku super-diagonals.
          interface gbmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,kl,ku,lda,m,n
                    character,intent(in) :: trans
                    complex(sp),intent(inout) :: y(*)
               end subroutine cgbmv
#else
               module procedure la_cgbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,kl,ku,lda,m,n
                    character,intent(in) :: trans
                    real(dp),intent(inout) :: y(*)
               end subroutine dgbmv
#else
               module procedure la_dgbmv
#endif
#ifdef LA_WITH_QP
               module procedure la_qgbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,kl,ku,lda,m,n
                    character,intent(in) :: trans
                    real(sp),intent(inout) :: y(*)
               end subroutine sgbmv
#else
               module procedure la_sgbmv
#endif
#ifdef LA_WITH_QP
               module procedure la_wgbmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xgbmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ygbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,kl,ku,lda,m,n
                    character,intent(in) :: trans
                    complex(dp),intent(inout) :: y(*)
               end subroutine zgbmv
#else
               module procedure la_zgbmv
#endif
          end interface gbmv

          !> GEMM:  performs one of the matrix-matrix operations
          !> C := alpha*op( A )*op( B ) + beta*C,
          !> where  op( X ) is one of
          !> op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
          !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
          !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
          interface gemm
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
                    character,intent(in) :: transa,transb
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine cgemm
#else
               module procedure la_cgemm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
                    character,intent(in) :: transa,transb
                    real(dp),intent(inout) :: c(ldc,*)
               end subroutine dgemm
#else
               module procedure la_dgemm
#endif
#ifdef LA_WITH_QP
               module procedure la_qgemm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
                    character,intent(in) :: transa,transb
                    real(sp),intent(inout) :: c(ldc,*)
               end subroutine sgemm
#else
               module procedure la_sgemm
#endif
#ifdef LA_WITH_QP
               module procedure la_wgemm
#endif
#ifdef LA_WITH_XDP
               module procedure la_xgemm
#endif
#ifdef LA_WITH_XDP
               module procedure la_ygemm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
                    character,intent(in) :: transa,transb
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zgemm
#else
               module procedure la_zgemm
#endif
          end interface gemm

          !> GEMV: performs one of the matrix-vector operations
          !> y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
          !> y := alpha*A**H*x + beta*y,
          !> where alpha and beta are scalars, x and y are vectors and A is an
          !> m by n matrix.
          interface gemv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    character,intent(in) :: trans
                    complex(sp),intent(inout) :: y(*)
               end subroutine cgemv
#else
               module procedure la_cgemv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    character,intent(in) :: trans
                    real(dp),intent(inout) :: y(*)
               end subroutine dgemv
#else
               module procedure la_dgemv
#endif
#ifdef LA_WITH_QP
               module procedure la_qgemv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    character,intent(in) :: trans
                    real(sp),intent(inout) :: y(*)
               end subroutine sgemv
#else
               module procedure la_sgemv
#endif
#ifdef LA_WITH_QP
               module procedure la_wgemv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xgemv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ygemv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    character,intent(in) :: trans
                    complex(dp),intent(inout) :: y(*)
               end subroutine zgemv
#else
               module procedure la_zgemv
#endif
          end interface gemv

          !> GER:   performs the rank 1 operation
          !> A := alpha*x*y**T + A,
          !> where alpha is a scalar, x is an m element vector, y is an n element
          !> vector and A is an m by n matrix.
          interface ger
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    real(dp),intent(inout) :: a(lda,*)
               end subroutine dger
#else
               module procedure la_dger
#endif
#ifdef LA_WITH_QP
               module procedure la_qger
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sger(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    real(sp),intent(inout) :: a(lda,*)
               end subroutine sger
#else
               module procedure la_sger
#endif
#ifdef LA_WITH_XDP
               module procedure la_xger
#endif
          end interface ger

          !> GERC:  performs the rank 1 operation
          !> A := alpha*x*y**H + A,
          !> where alpha is a scalar, x is an m element vector, y is an n element
          !> vector and A is an m by n matrix.
          interface gerc
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cgerc(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    complex(sp),intent(inout) :: a(lda,*)
               end subroutine cgerc
#else
               module procedure la_cgerc
#endif
#ifdef LA_WITH_QP
               module procedure la_wgerc
#endif
#ifdef LA_WITH_XDP
               module procedure la_ygerc
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zgerc(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    complex(dp),intent(inout) :: a(lda,*)
               end subroutine zgerc
#else
               module procedure la_zgerc
#endif
          end interface gerc

          !> GERU:  performs the rank 1 operation
          !> A := alpha*x*y**T + A,
          !> where alpha is a scalar, x is an m element vector, y is an n element
          !> vector and A is an m by n matrix.
          interface geru
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cgeru(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    complex(sp),intent(inout) :: a(lda,*)
               end subroutine cgeru
#else
               module procedure la_cgeru
#endif
#ifdef LA_WITH_QP
               module procedure la_wgeru
#endif
#ifdef LA_WITH_XDP
               module procedure la_ygeru
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zgeru(m,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,m,n
                    complex(dp),intent(inout) :: a(lda,*)
               end subroutine zgeru
#else
               module procedure la_zgeru
#endif
          end interface geru

          !> HBMV:  performs the matrix-vector  operation
          !> y := alpha*A*x + beta*y,
          !> where alpha and beta are scalars, x and y are n element vectors and
          !> A is an n by n hermitian band matrix, with k super-diagonals.
          interface hbmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine chbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,k,lda,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: y(*)
               end subroutine chbmv
#else
               module procedure la_chbmv
#endif
#ifdef LA_WITH_QP
               module procedure la_whbmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_yhbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zhbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,k,lda,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: y(*)
               end subroutine zhbmv
#else
               module procedure la_zhbmv
#endif
          end interface hbmv

          !> HEMM:  performs one of the matrix-matrix operations
          !> C := alpha*A*B + beta*C,
          !> or
          !> C := alpha*B*A + beta*C,
          !> where alpha and beta are scalars, A is an hermitian matrix and  B and
          !> C are m by n matrices.
          interface hemm
#ifdef LA_EXTERNAL_BLAS
               pure subroutine chemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: lda,ldb,ldc,m,n
                    character,intent(in) :: side,uplo
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine chemm
#else
               module procedure la_chemm
#endif
#ifdef LA_WITH_QP
               module procedure la_whemm
#endif
#ifdef LA_WITH_XDP
               module procedure la_yhemm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: lda,ldb,ldc,m,n
                    character,intent(in) :: side,uplo
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zhemm
#else
               module procedure la_zhemm
#endif
          end interface hemm

          !> HEMV:  performs the matrix-vector  operation
          !> y := alpha*A*x + beta*y,
          !> where alpha and beta are scalars, x and y are n element vectors and
          !> A is an n by n hermitian matrix.
          interface hemv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine chemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: y(*)
               end subroutine chemv
#else
               module procedure la_chemv
#endif
#ifdef LA_WITH_QP
               module procedure la_whemv
#endif
#ifdef LA_WITH_XDP
               module procedure la_yhemv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: y(*)
               end subroutine zhemv
#else
               module procedure la_zhemv
#endif
          end interface hemv

          !> HER:   performs the hermitian rank 1 operation
          !> A := alpha*x*x**H + A,
          !> where alpha is a real scalar, x is an n element vector and A is an
          !> n by n hermitian matrix.
          interface her
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cher(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: a(lda,*)
                    complex(sp),intent(in) :: x(*)
               end subroutine cher
#else
               module procedure la_cher
#endif
#ifdef LA_WITH_QP
               module procedure la_wher
#endif
#ifdef LA_WITH_XDP
               module procedure la_yher
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zher(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: a(lda,*)
                    complex(dp),intent(in) :: x(*)
               end subroutine zher
#else
               module procedure la_zher
#endif
          end interface her

          !> HER2:  performs the hermitian rank 2 operation
          !> A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
          !> where alpha is a scalar, x and y are n element vectors and A is an n
          !> by n hermitian matrix.
          interface her2
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cher2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: a(lda,*)
               end subroutine cher2
#else
               module procedure la_cher2
#endif
#ifdef LA_WITH_QP
               module procedure la_wher2
#endif
#ifdef LA_WITH_XDP
               module procedure la_yher2
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: a(lda,*)
               end subroutine zher2
#else
               module procedure la_zher2
#endif
          end interface her2

          !> HER2K:  performs one of the hermitian rank 2k operations
          !> C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
          !> or
          !> C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
          !> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
          !> hermitian matrix and  A and B  are  n by k matrices in the first case
          !> and  k by n  matrices in the second case.
          interface her2k
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,a(lda,*),b(ldb,*)
                    real(sp),intent(in) :: beta
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine cher2k
#else
               module procedure la_cher2k
#endif
#ifdef LA_WITH_QP
               module procedure la_wher2k
#endif
#ifdef LA_WITH_XDP
               module procedure la_yher2k
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,a(lda,*),b(ldb,*)
                    real(dp),intent(in) :: beta
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zher2k
#else
               module procedure la_zher2k
#endif
          end interface her2k

          !> HERK:  performs one of the hermitian rank k operations
          !> C := alpha*A*A**H + beta*C,
          !> or
          !> C := alpha*A**H*A + beta*C,
          !> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
          !> matrix and  A  is an  n by k  matrix in the  first case and a  k by n
          !> matrix in the second case.
          interface herk
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta
                    integer(ilp),intent(in) :: k,lda,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(sp),intent(in) :: a(lda,*)
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine cherk
#else
               module procedure la_cherk
#endif
#ifdef LA_WITH_QP
               module procedure la_wherk
#endif
#ifdef LA_WITH_XDP
               module procedure la_yherk
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta
                    integer(ilp),intent(in) :: k,lda,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(dp),intent(in) :: a(lda,*)
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zherk
#else
               module procedure la_zherk
#endif
          end interface herk

          !> HPMV:  performs the matrix-vector operation
          !> y := alpha*A*x + beta*y,
          !> where alpha and beta are scalars, x and y are n element vectors and
          !> A is an n by n hermitian matrix, supplied in packed form.
          interface hpmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine chpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,ap(*),x(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: y(*)
               end subroutine chpmv
#else
               module procedure la_chpmv
#endif
#ifdef LA_WITH_QP
               module procedure la_whpmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_yhpmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,ap(*),x(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: y(*)
               end subroutine zhpmv
#else
               module procedure la_zhpmv
#endif
          end interface hpmv

          !> HPR:    performs the hermitian rank 1 operation
          !> A := alpha*x*x**H + A,
          !> where alpha is a real scalar, x is an n element vector and A is an
          !> n by n hermitian matrix, supplied in packed form.
          interface hpr
#ifdef LA_EXTERNAL_BLAS
               pure subroutine chpr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: ap(*)
                    complex(sp),intent(in) :: x(*)
               end subroutine chpr
#else
               module procedure la_chpr
#endif
#ifdef LA_WITH_QP
               module procedure la_whpr
#endif
#ifdef LA_WITH_XDP
               module procedure la_yhpr
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zhpr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: ap(*)
                    complex(dp),intent(in) :: x(*)
               end subroutine zhpr
#else
               module procedure la_zhpr
#endif
          end interface hpr

          !> HPR2:  performs the hermitian rank 2 operation
          !> A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
          !> where alpha is a scalar, x and y are n element vectors and A is an
          !> n by n hermitian matrix, supplied in packed form.
          interface hpr2
#ifdef LA_EXTERNAL_BLAS
               pure subroutine chpr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    complex(sp),intent(inout) :: ap(*)
               end subroutine chpr2
#else
               module procedure la_chpr2
#endif
#ifdef LA_WITH_QP
               module procedure la_whpr2
#endif
#ifdef LA_WITH_XDP
               module procedure la_yhpr2
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    complex(dp),intent(inout) :: ap(*)
               end subroutine zhpr2
#else
               module procedure la_zhpr2
#endif
          end interface hpr2

          !> !
          !>
          !> NRM2: returns the euclidean norm of a vector via the function
          !> name, so that
          !> NRM2 := sqrt( x'*x )
          interface nrm2
#ifdef LA_EXTERNAL_BLAS
               pure real(dp) function dnrm2(n,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    real(dp),intent(in) :: x(*)
               end function dnrm2
#else
               module procedure la_dnrm2
#endif
#ifdef LA_WITH_QP
               module procedure la_qnrm2
#endif
#ifdef LA_EXTERNAL_BLAS
               pure real(sp) function snrm2(n,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    real(sp),intent(in) :: x(*)
               end function snrm2
#else
               module procedure la_snrm2
#endif
#ifdef LA_WITH_XDP
               module procedure la_xnrm2
#endif
          end interface nrm2

          !> ROT: applies a plane rotation.
          interface rot
#ifdef LA_EXTERNAL_BLAS
               pure subroutine drot(n,dx,incx,dy,incy,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: c,s
                    integer(ilp),intent(in) :: incx,incy,n
                    real(dp),intent(inout) :: dx(*),dy(*)
               end subroutine drot
#else
               module procedure la_drot
#endif
#ifdef LA_WITH_QP
               module procedure la_qrot
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine srot(n,sx,incx,sy,incy,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: c,s
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(inout) :: sx(*),sy(*)
               end subroutine srot
#else
               module procedure la_srot
#endif
#ifdef LA_WITH_XDP
               module procedure la_xrot
#endif
          end interface rot

          !> !
          !>
          !> The computation uses the formulas
          !> |x| = sqrt( Re(x)**2 + Im(x)**2 )
          !> sgn(x) = x / |x|  if x /= 0
          !> = 1        if x  = 0
          !> c = |a| / sqrt(|a|**2 + |b|**2)
          !> s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
          !> When a and b are real and r /= 0, the formulas simplify to
          !> r = sgn(a)*sqrt(|a|**2 + |b|**2)
          !> c = a / r
          !> s = b / r
          !> the same as in SROTG when |a| > |b|.  When |b| >= |a|, the
          !> sign of c and s will be different from those computed by SROTG
          !> if the signs of a and b are not the same.
          interface rotg
#ifdef LA_EXTERNAL_BLAS
               pure subroutine crotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(out) :: c
                    complex(sp),intent(inout) :: a
                    complex(sp),intent(in) :: b
                    complex(sp),intent(out) :: s
               end subroutine crotg
#else
               module procedure la_crotg
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine drotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(inout) :: a,b
                    real(dp),intent(out) :: c,s
               end subroutine drotg
#else
               module procedure la_drotg
#endif
#ifdef LA_WITH_QP
               module procedure la_qrotg
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine srotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(inout) :: a,b
                    real(sp),intent(out) :: c,s
               end subroutine srotg
#else
               module procedure la_srotg
#endif
#ifdef LA_WITH_QP
               module procedure la_wrotg
#endif
#ifdef LA_WITH_XDP
               module procedure la_xrotg
#endif
#ifdef LA_WITH_XDP
               module procedure la_yrotg
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zrotg(a,b,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(out) :: c
                    complex(dp),intent(inout) :: a
                    complex(dp),intent(in) :: b
                    complex(dp),intent(out) :: s
               end subroutine zrotg
#else
               module procedure la_zrotg
#endif
          end interface rotg

          !> APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
          !> (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
          !> (DY**T)
          !> DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
          !> LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
          !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
          !> DFLAG=-1._dp     DFLAG=0._dp        DFLAG=1._dp     DFLAG=-2.D0
          !> (DH11  DH12)    (1._dp  DH12)    (DH11  1._dp)    (1._dp  0._dp)
          !> H=(          )    (          )    (          )    (          )
          !> (DH21  DH22),   (DH21  1._dp),   (-1._dp DH22),   (0._dp  1._dp).
          !> SEE ROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
          interface rotm
#ifdef LA_EXTERNAL_BLAS
               pure subroutine drotm(n,dx,incx,dy,incy,dparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(dp),intent(in) :: dparam(5)
                    real(dp),intent(inout) :: dx(*),dy(*)
               end subroutine drotm
#else
               module procedure la_drotm
#endif
#ifdef LA_WITH_QP
               module procedure la_qrotm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine srotm(n,sx,incx,sy,incy,sparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(in) :: sparam(5)
                    real(sp),intent(inout) :: sx(*),sy(*)
               end subroutine srotm
#else
               module procedure la_srotm
#endif
#ifdef LA_WITH_XDP
               module procedure la_xrotm
#endif
          end interface rotm

          !> CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
          !> THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(DD1)*DX1,SQRT(DD2)    DY2)**T.
          !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
          !> DFLAG=-1._dp     DFLAG=0._dp        DFLAG=1._dp     DFLAG=-2.D0
          !> (DH11  DH12)    (1._dp  DH12)    (DH11  1._dp)    (1._dp  0._dp)
          !> H=(          )    (          )    (          )    (          )
          !> (DH21  DH22),   (DH21  1._dp),   (-1._dp DH22),   (0._dp  1._dp).
          !> LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
          !> RESPECTIVELY. (VALUES OF 1._dp, -1._dp, OR 0._dp IMPLIED BY THE
          !> VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
          !> THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
          !> INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
          !> OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
          interface rotmg
#ifdef LA_EXTERNAL_BLAS
               pure subroutine drotmg(dd1,dd2,dx1,dy1,dparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(inout) :: dd1,dd2,dx1
                    real(dp),intent(in) :: dy1
                    real(dp),intent(out) :: dparam(5)
               end subroutine drotmg
#else
               module procedure la_drotmg
#endif
#ifdef LA_WITH_QP
               module procedure la_qrotmg
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine srotmg(sd1,sd2,sx1,sy1,sparam)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(inout) :: sd1,sd2,sx1
                    real(sp),intent(in) :: sy1
                    real(sp),intent(out) :: sparam(5)
               end subroutine srotmg
#else
               module procedure la_srotmg
#endif
#ifdef LA_WITH_XDP
               module procedure la_xrotmg
#endif
          end interface rotmg

          !> SBMV:  performs the matrix-vector  operation
          !> y := alpha*A*x + beta*y,
          !> where alpha and beta are scalars, x and y are n element vectors and
          !> A is an n by n symmetric band matrix, with k super-diagonals.
          interface sbmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,k,lda,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: y(*)
               end subroutine dsbmv
#else
               module procedure la_dsbmv
#endif
#ifdef LA_WITH_QP
               module procedure la_qsbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,k,lda,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: y(*)
               end subroutine ssbmv
#else
               module procedure la_ssbmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsbmv
#endif
          end interface sbmv

          !> SCAL: scales a vector by a constant.
          interface scal
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cscal(n,ca,cx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: ca
                    integer(ilp),intent(in) :: incx,n
                    complex(sp),intent(inout) :: cx(*)
               end subroutine cscal
#else
               module procedure la_cscal
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dscal(n,da,dx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: da
                    integer(ilp),intent(in) :: incx,n
                    real(dp),intent(inout) :: dx(*)
               end subroutine dscal
#else
               module procedure la_dscal
#endif
#ifdef LA_WITH_QP
               module procedure la_qscal
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sscal(n,sa,sx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: sa
                    integer(ilp),intent(in) :: incx,n
                    real(sp),intent(inout) :: sx(*)
               end subroutine sscal
#else
               module procedure la_sscal
#endif
#ifdef LA_WITH_QP
               module procedure la_wscal
#endif
#ifdef LA_WITH_XDP
               module procedure la_xscal
#endif
#ifdef LA_WITH_XDP
               module procedure la_yscal
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zscal(n,za,zx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: za
                    integer(ilp),intent(in) :: incx,n
                    complex(dp),intent(inout) :: zx(*)
               end subroutine zscal
#else
               module procedure la_zscal
#endif
          end interface scal

          !> Compute the inner product of two vectors with extended
          !> precision accumulation and result.
          !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
          !> SDOT: = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
          !> where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
          !> defined in a similar way using INCY.
          interface sdot
#ifdef LA_EXTERNAL_BLAS
               pure real(dp) function dsdot(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(in) :: sx(*),sy(*)
               end function dsdot
#else
               module procedure la_dsdot
#endif
#ifdef LA_WITH_QP
               module procedure la_qddot
#endif
#ifdef LA_WITH_XDP
               module procedure la_xddot
#endif
          end interface sdot

          !> SPMV:  performs the matrix-vector operation
          !> y := alpha*A*x + beta*y,
          !> where alpha and beta are scalars, x and y are n element vectors and
          !> A is an n by n symmetric matrix, supplied in packed form.
          interface spmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,ap(*),x(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: y(*)
               end subroutine dspmv
#else
               module procedure la_dspmv
#endif
#ifdef LA_WITH_QP
               module procedure la_qspmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,ap(*),x(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: y(*)
               end subroutine sspmv
#else
               module procedure la_sspmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xspmv
#endif
          end interface spmv

          !> SPR:    performs the symmetric rank 1 operation
          !> A := alpha*x*x**T + A,
          !> where alpha is a real scalar, x is an n element vector and A is an
          !> n by n symmetric matrix, supplied in packed form.
          interface spr
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dspr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,x(*)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: ap(*)
               end subroutine dspr
#else
               module procedure la_dspr
#endif
#ifdef LA_WITH_QP
               module procedure la_qspr
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sspr(uplo,n,alpha,x,incx,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,x(*)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: ap(*)
               end subroutine sspr
#else
               module procedure la_sspr
#endif
#ifdef LA_WITH_XDP
               module procedure la_xspr
#endif
          end interface spr

          !> SPR2:  performs the symmetric rank 2 operation
          !> A := alpha*x*y**T + alpha*y*x**T + A,
          !> where alpha is a scalar, x and y are n element vectors and A is an
          !> n by n symmetric matrix, supplied in packed form.
          interface spr2
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dspr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: ap(*)
               end subroutine dspr2
#else
               module procedure la_dspr2
#endif
#ifdef LA_WITH_QP
               module procedure la_qspr2
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sspr2(uplo,n,alpha,x,incx,y,incy,ap)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: ap(*)
               end subroutine sspr2
#else
               module procedure la_sspr2
#endif
#ifdef LA_WITH_XDP
               module procedure la_xspr2
#endif
          end interface spr2

          !> SROT: applies a plane rotation, where the cos and sin (c and s) are real
          !> and the vectors cx and cy are complex.
          !> jack dongarra, linpack, 3/11/78.
          interface srot
#ifdef LA_EXTERNAL_BLAS
               pure subroutine csrot(n,cx,incx,cy,incy,c,s)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(in) :: c,s
                    complex(sp),intent(inout) :: cx(*),cy(*)
               end subroutine csrot
#else
               module procedure la_csrot
#endif
          end interface srot

          !> SSCAL: scales a complex vector by a real constant.
          interface sscal
#ifdef LA_EXTERNAL_BLAS
               pure subroutine csscal(n,sa,cx,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: sa
                    integer(ilp),intent(in) :: incx,n
                    complex(sp),intent(inout) :: cx(*)
               end subroutine csscal
#else
               module procedure la_csscal
#endif
          end interface sscal

          !> SWAP: interchanges two vectors.
          interface swap
#ifdef LA_EXTERNAL_BLAS
               pure subroutine cswap(n,cx,incx,cy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(sp),intent(inout) :: cx(*),cy(*)
               end subroutine cswap
#else
               module procedure la_cswap
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dswap(n,dx,incx,dy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(dp),intent(inout) :: dx(*),dy(*)
               end subroutine dswap
#else
               module procedure la_dswap
#endif
#ifdef LA_WITH_QP
               module procedure la_qswap
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine sswap(n,sx,incx,sy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    real(sp),intent(inout) :: sx(*),sy(*)
               end subroutine sswap
#else
               module procedure la_sswap
#endif
#ifdef LA_WITH_QP
               module procedure la_wswap
#endif
#ifdef LA_WITH_XDP
               module procedure la_xswap
#endif
#ifdef LA_WITH_XDP
               module procedure la_yswap
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zswap(n,zx,incx,zy,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,incy,n
                    complex(dp),intent(inout) :: zx(*),zy(*)
               end subroutine zswap
#else
               module procedure la_zswap
#endif
          end interface swap

          !> SYMM:  performs one of the matrix-matrix operations
          !> C := alpha*A*B + beta*C,
          !> or
          !> C := alpha*B*A + beta*C,
          !> where  alpha and beta are scalars, A is a symmetric matrix and  B and
          !> C are m by n matrices.
          interface symm
#ifdef LA_EXTERNAL_BLAS
               pure subroutine csymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: lda,ldb,ldc,m,n
                    character,intent(in) :: side,uplo
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine csymm
#else
               module procedure la_csymm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: lda,ldb,ldc,m,n
                    character,intent(in) :: side,uplo
                    real(dp),intent(inout) :: c(ldc,*)
               end subroutine dsymm
#else
               module procedure la_dsymm
#endif
#ifdef LA_WITH_QP
               module procedure la_qsymm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: lda,ldb,ldc,m,n
                    character,intent(in) :: side,uplo
                    real(sp),intent(inout) :: c(ldc,*)
               end subroutine ssymm
#else
               module procedure la_ssymm
#endif
#ifdef LA_WITH_QP
               module procedure la_wsymm
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsymm
#endif
#ifdef LA_WITH_XDP
               module procedure la_ysymm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: lda,ldb,ldc,m,n
                    character,intent(in) :: side,uplo
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zsymm
#else
               module procedure la_zsymm
#endif
          end interface symm

          !> SYMV:  performs the matrix-vector  operation
          !> y := alpha*A*x + beta*y,
          !> where alpha and beta are scalars, x and y are n element vectors and
          !> A is an n by n symmetric matrix.
          interface symv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: y(*)
               end subroutine dsymv
#else
               module procedure la_dsymv
#endif
#ifdef LA_WITH_QP
               module procedure la_qsymv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),x(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: y(*)
               end subroutine ssymv
#else
               module procedure la_ssymv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsymv
#endif
          end interface symv

          !> SYR:   performs the symmetric rank 1 operation
          !> A := alpha*x*x**T + A,
          !> where alpha is a real scalar, x is an n element vector and A is an
          !> n by n symmetric matrix.
          interface syr
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsyr(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,x(*)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: a(lda,*)
               end subroutine dsyr
#else
               module procedure la_dsyr
#endif
#ifdef LA_WITH_QP
               module procedure la_qsyr
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssyr(uplo,n,alpha,x,incx,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,x(*)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: a(lda,*)
               end subroutine ssyr
#else
               module procedure la_ssyr
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsyr
#endif
          end interface syr

          !> SYR2:  performs the symmetric rank 2 operation
          !> A := alpha*x*y**T + alpha*y*x**T + A,
          !> where alpha is a scalar, x and y are n element vectors and A is an n
          !> by n symmetric matrix.
          interface syr2
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    real(dp),intent(inout) :: a(lda,*)
               end subroutine dsyr2
#else
               module procedure la_dsyr2
#endif
#ifdef LA_WITH_QP
               module procedure la_qsyr2
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,x(*),y(*)
                    integer(ilp),intent(in) :: incx,incy,lda,n
                    character,intent(in) :: uplo
                    real(sp),intent(inout) :: a(lda,*)
               end subroutine ssyr2
#else
               module procedure la_ssyr2
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsyr2
#endif
          end interface syr2

          !> SYR2K:  performs one of the symmetric rank 2k operations
          !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
          !> or
          !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
          !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
          !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
          !> matrices in the second case.
          interface syr2k
#ifdef LA_EXTERNAL_BLAS
               pure subroutine csyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine csyr2k
#else
               module procedure la_csyr2k
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,n
                    character,intent(in) :: trans,uplo
                    real(dp),intent(inout) :: c(ldc,*)
               end subroutine dsyr2k
#else
               module procedure la_dsyr2k
#endif
#ifdef LA_WITH_QP
               module procedure la_qsyr2k
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,n
                    character,intent(in) :: trans,uplo
                    real(sp),intent(inout) :: c(ldc,*)
               end subroutine ssyr2k
#else
               module procedure la_ssyr2k
#endif
#ifdef LA_WITH_QP
               module procedure la_wsyr2k
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsyr2k
#endif
#ifdef LA_WITH_XDP
               module procedure la_ysyr2k
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
                    integer(ilp),intent(in) :: k,lda,ldb,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zsyr2k
#else
               module procedure la_zsyr2k
#endif
          end interface syr2k

          !> SYRK:  performs one of the symmetric rank k operations
          !> C := alpha*A*A**T + beta*C,
          !> or
          !> C := alpha*A**T*A + beta*C,
          !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
          !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
          !> in the second case.
          interface syrk
#ifdef LA_EXTERNAL_BLAS
               pure subroutine csyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,beta,a(lda,*)
                    integer(ilp),intent(in) :: k,lda,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(sp),intent(inout) :: c(ldc,*)
               end subroutine csyrk
#else
               module procedure la_csyrk
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,beta,a(lda,*)
                    integer(ilp),intent(in) :: k,lda,ldc,n
                    character,intent(in) :: trans,uplo
                    real(dp),intent(inout) :: c(ldc,*)
               end subroutine dsyrk
#else
               module procedure la_dsyrk
#endif
#ifdef LA_WITH_QP
               module procedure la_qsyrk
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ssyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,beta,a(lda,*)
                    integer(ilp),intent(in) :: k,lda,ldc,n
                    character,intent(in) :: trans,uplo
                    real(sp),intent(inout) :: c(ldc,*)
               end subroutine ssyrk
#else
               module procedure la_ssyrk
#endif
#ifdef LA_WITH_QP
               module procedure la_wsyrk
#endif
#ifdef LA_WITH_XDP
               module procedure la_xsyrk
#endif
#ifdef LA_WITH_XDP
               module procedure la_ysyrk
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,beta,a(lda,*)
                    integer(ilp),intent(in) :: k,lda,ldc,n
                    character,intent(in) :: trans,uplo
                    complex(dp),intent(inout) :: c(ldc,*)
               end subroutine zsyrk
#else
               module procedure la_zsyrk
#endif
          end interface syrk

          !> TBMV:  performs one of the matrix-vector operations
          !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
          !> where x is an n element vector and  A is an n by n unit, or non-unit,
          !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.
          interface tbmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(sp),intent(in) :: a(lda,*)
                    complex(sp),intent(inout) :: x(*)
               end subroutine ctbmv
#else
               module procedure la_ctbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(dp),intent(in) :: a(lda,*)
                    real(dp),intent(inout) :: x(*)
               end subroutine dtbmv
#else
               module procedure la_dtbmv
#endif
#ifdef LA_WITH_QP
               module procedure la_qtbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine stbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(sp),intent(in) :: a(lda,*)
                    real(sp),intent(inout) :: x(*)
               end subroutine stbmv
#else
               module procedure la_stbmv
#endif
#ifdef LA_WITH_QP
               module procedure la_wtbmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtbmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytbmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztbmv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(dp),intent(in) :: a(lda,*)
                    complex(dp),intent(inout) :: x(*)
               end subroutine ztbmv
#else
               module procedure la_ztbmv
#endif
          end interface tbmv

          !> TBSV:  solves one of the systems of equations
          !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
          !> where b and x are n element vectors and A is an n by n unit, or
          !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
          !> diagonals.
          !> No test for singularity or near-singularity is included in this
          !> routine. Such tests must be performed before calling this routine.
          interface tbsv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(sp),intent(in) :: a(lda,*)
                    complex(sp),intent(inout) :: x(*)
               end subroutine ctbsv
#else
               module procedure la_ctbsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(dp),intent(in) :: a(lda,*)
                    real(dp),intent(inout) :: x(*)
               end subroutine dtbsv
#else
               module procedure la_dtbsv
#endif
#ifdef LA_WITH_QP
               module procedure la_qtbsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine stbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(sp),intent(in) :: a(lda,*)
                    real(sp),intent(inout) :: x(*)
               end subroutine stbsv
#else
               module procedure la_stbsv
#endif
#ifdef LA_WITH_QP
               module procedure la_wtbsv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtbsv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytbsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztbsv(uplo,trans,diag,n,k,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,k,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(dp),intent(in) :: a(lda,*)
                    complex(dp),intent(inout) :: x(*)
               end subroutine ztbsv
#else
               module procedure la_ztbsv
#endif
          end interface tbsv

          !> TPMV:  performs one of the matrix-vector operations
          !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
          !> where x is an n element vector and  A is an n by n unit, or non-unit,
          !> upper or lower triangular matrix, supplied in packed form.
          interface tpmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    complex(sp),intent(in) :: ap(*)
                    complex(sp),intent(inout) :: x(*)
               end subroutine ctpmv
#else
               module procedure la_ctpmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    real(dp),intent(in) :: ap(*)
                    real(dp),intent(inout) :: x(*)
               end subroutine dtpmv
#else
               module procedure la_dtpmv
#endif
#ifdef LA_WITH_QP
               module procedure la_qtpmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine stpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    real(sp),intent(in) :: ap(*)
                    real(sp),intent(inout) :: x(*)
               end subroutine stpmv
#else
               module procedure la_stpmv
#endif
#ifdef LA_WITH_QP
               module procedure la_wtpmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtpmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytpmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztpmv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    complex(dp),intent(in) :: ap(*)
                    complex(dp),intent(inout) :: x(*)
               end subroutine ztpmv
#else
               module procedure la_ztpmv
#endif
          end interface tpmv

          !> TPSV:  solves one of the systems of equations
          !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
          !> where b and x are n element vectors and A is an n by n unit, or
          !> non-unit, upper or lower triangular matrix, supplied in packed form.
          !> No test for singularity or near-singularity is included in this
          !> routine. Such tests must be performed before calling this routine.
          interface tpsv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    complex(sp),intent(in) :: ap(*)
                    complex(sp),intent(inout) :: x(*)
               end subroutine ctpsv
#else
               module procedure la_ctpsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    real(dp),intent(in) :: ap(*)
                    real(dp),intent(inout) :: x(*)
               end subroutine dtpsv
#else
               module procedure la_dtpsv
#endif
#ifdef LA_WITH_QP
               module procedure la_qtpsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine stpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    real(sp),intent(in) :: ap(*)
                    real(sp),intent(inout) :: x(*)
               end subroutine stpsv
#else
               module procedure la_stpsv
#endif
#ifdef LA_WITH_QP
               module procedure la_wtpsv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtpsv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytpsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztpsv(uplo,trans,diag,n,ap,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,n
                    character,intent(in) :: diag,trans,uplo
                    complex(dp),intent(in) :: ap(*)
                    complex(dp),intent(inout) :: x(*)
               end subroutine ztpsv
#else
               module procedure la_ztpsv
#endif
          end interface tpsv

          !> TRMM:  performs one of the matrix-matrix operations
          !> B := alpha*op( A )*B,   or   B := alpha*B*op( A )
          !> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
          !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
          !> op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
          interface trmm
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    complex(sp),intent(inout) :: b(ldb,*)
               end subroutine ctrmm
#else
               module procedure la_ctrmm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    real(dp),intent(inout) :: b(ldb,*)
               end subroutine dtrmm
#else
               module procedure la_dtrmm
#endif
#ifdef LA_WITH_QP
               module procedure la_qtrmm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine strmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    real(sp),intent(inout) :: b(ldb,*)
               end subroutine strmm
#else
               module procedure la_strmm
#endif
#ifdef LA_WITH_QP
               module procedure la_wtrmm
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtrmm
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytrmm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    complex(dp),intent(inout) :: b(ldb,*)
               end subroutine ztrmm
#else
               module procedure la_ztrmm
#endif
          end interface trmm

          !> TRMV:  performs one of the matrix-vector operations
          !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
          !> where x is an n element vector and  A is an n by n unit, or non-unit,
          !> upper or lower triangular matrix.
          interface trmv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctrmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(sp),intent(in) :: a(lda,*)
                    complex(sp),intent(inout) :: x(*)
               end subroutine ctrmv
#else
               module procedure la_ctrmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtrmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(dp),intent(in) :: a(lda,*)
                    real(dp),intent(inout) :: x(*)
               end subroutine dtrmv
#else
               module procedure la_dtrmv
#endif
#ifdef LA_WITH_QP
               module procedure la_qtrmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine strmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(sp),intent(in) :: a(lda,*)
                    real(sp),intent(inout) :: x(*)
               end subroutine strmv
#else
               module procedure la_strmv
#endif
#ifdef LA_WITH_QP
               module procedure la_wtrmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtrmv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytrmv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztrmv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(dp),intent(in) :: a(lda,*)
                    complex(dp),intent(inout) :: x(*)
               end subroutine ztrmv
#else
               module procedure la_ztrmv
#endif
          end interface trmv

          !> TRSM:  solves one of the matrix equations
          !> op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
          !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
          !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
          !> op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
          !> The matrix X is overwritten on B.
          interface trsm
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(sp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    complex(sp),intent(inout) :: b(ldb,*)
               end subroutine ctrsm
#else
               module procedure la_ctrsm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(dp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    real(dp),intent(inout) :: b(ldb,*)
               end subroutine dtrsm
#else
               module procedure la_dtrsm
#endif
#ifdef LA_WITH_QP
               module procedure la_qtrsm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine strsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    real(sp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    real(sp),intent(inout) :: b(ldb,*)
               end subroutine strsm
#else
               module procedure la_strsm
#endif
#ifdef LA_WITH_QP
               module procedure la_wtrsm
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtrsm
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytrsm
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    complex(dp),intent(in) :: alpha,a(lda,*)
                    integer(ilp),intent(in) :: lda,ldb,m,n
                    character,intent(in) :: diag,side,transa,uplo
                    complex(dp),intent(inout) :: b(ldb,*)
               end subroutine ztrsm
#else
               module procedure la_ztrsm
#endif
          end interface trsm

          !> TRSV:  solves one of the systems of equations
          !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
          !> where b and x are n element vectors and A is an n by n unit, or
          !> non-unit, upper or lower triangular matrix.
          !> No test for singularity or near-singularity is included in this
          !> routine. Such tests must be performed before calling this routine.
          interface trsv
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ctrsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(sp),intent(in) :: a(lda,*)
                    complex(sp),intent(inout) :: x(*)
               end subroutine ctrsv
#else
               module procedure la_ctrsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(dp),intent(in) :: a(lda,*)
                    real(dp),intent(inout) :: x(*)
               end subroutine dtrsv
#else
               module procedure la_dtrsv
#endif
#ifdef LA_WITH_QP
               module procedure la_qtrsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine strsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    real(sp),intent(in) :: a(lda,*)
                    real(sp),intent(inout) :: x(*)
               end subroutine strsv
#else
               module procedure la_strsv
#endif
#ifdef LA_WITH_QP
               module procedure la_wtrsv
#endif
#ifdef LA_WITH_XDP
               module procedure la_xtrsv
#endif
#ifdef LA_WITH_XDP
               module procedure la_ytrsv
#endif
#ifdef LA_EXTERNAL_BLAS
               pure subroutine ztrsv(uplo,trans,diag,n,a,lda,x,incx)
                    import sp,dp,qp,ilp,lk
                    implicit none(type,external)
                    integer(ilp),intent(in) :: incx,lda,n
                    character,intent(in) :: diag,trans,uplo
                    complex(dp),intent(in) :: a(lda,*)
                    complex(dp),intent(inout) :: x(*)
               end subroutine ztrsv
#else
               module procedure la_ztrsv
#endif
          end interface trsv

end module la_blas
