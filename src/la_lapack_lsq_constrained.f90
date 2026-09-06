!> Constrained least squares: equality constraints and the general Gauss-Markov model
module la_lapack_lsq_constrained
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_tri
     use la_lapack_aux
     use la_lapack_orthogonal_factors_qr
     use la_lapack_solve_tri_comp
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sggglm
     public :: la_sgglse
     public :: la_dggglm
     public :: la_dgglse
     public :: la_qggglm
     public :: la_qgglse
     public :: la_cggglm
     public :: la_cgglse
     public :: la_zggglm
     public :: la_zgglse
     public :: la_wggglm
     public :: la_wgglse

     contains

     !> SGGGLM: solves a general Gauss-Markov linear model (GLM) problem:
     !> minimize || y ||_2   subject to   d = A*x + B*y
     !> x
     !> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     !> given N-vector. It is assumed that M <= N <= M+P, and
     !> rank(A) = M    and    rank( A B ) = N.
     !> Under these assumptions, the constrained equation is always
     !> consistent, and there is a unique solution x and a minimal 2-norm
     !> solution y, which is obtained using a generalized QR factorization
     !> of the matrices (A, B) given by
     !> A = Q*(R),   B = Q*T*Z.
     !> (0)
     !> In particular, if matrix B is square nonsingular, then the problem
     !> GLM is equivalent to the following weighted linear least squares
     !> problem
     !> minimize || inv(B)*(d-A*x) ||_2
     !> x
     !> where inv(B) denotes the inverse of B.

     pure subroutine la_sggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*),d(*)
           real(sp),intent(out) :: work(*),x(*),y(*)
        ! ===================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,lopt,lwkmin,lwkopt,nb,nb1,nb2,nb3,nb4,np
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           np = min(n,p)
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -2
           else if (p < 0 .or. p < n - m) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'SGEQRF',' ',n,m,-1,-1)
                 nb2 = la_ilaenv(1,'SGERQF',' ',n,m,-1,-1)
                 nb3 = la_ilaenv(1,'SORMQR',' ',n,m,p,-1)
                 nb4 = la_ilaenv(1,'SORMRQ',' ',n,m,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = m + np + max(n,p)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGGGLM',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              do i = 1,m
                 x(i) = zero
              end do
              do i = 1,p
                 y(i) = zero
              end do
              return
           end if
           ! compute the gqr factorization of matrices a and b:
                ! q**t*a = ( r11 ) m,    q**t*b*z**t = ( t11   t12 ) m
                         ! (  0  ) n-m                 (  0    t22 ) n-m
                            ! m                         m+p-n  n-m
           ! where r11 and t22 are upper triangular, and q and z are
           ! orthogonal.
           call la_sggqrf(n,m,p,a,lda,work,b,ldb,work(m + 1),work(m + np + 1),lwork - m - &
                     np,info)
           lopt = work(m + np + 1)
           ! update left-hand-side vector d = q**t*d = ( d1 ) m
                                                     ! ( d2 ) n-m
           call la_sormqr('LEFT','TRANSPOSE',n,1,m,a,lda,work,d,max(1,n),work(m + &
                     np + 1),lwork - m - np,info)
           lopt = max(lopt,int(work(m + np + 1),KIND=ilp))
           ! solve t22*y2 = d2 for y2
           if (n > m) then
              call la_strtrs('UPPER','NO TRANSPOSE','NON UNIT',n - m,1,b(m + 1,m + p - n + 1), &
                        ldb,d(m + 1),n - m,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              call la_scopy(n - m,d(m + 1),1,y(m + p - n + 1),1)
           end if
           ! set y1 = 0
           do i = 1,m + p - n
              y(i) = zero
           end do
           ! update d1 = d1 - t12*y2
           call la_sgemv('NO TRANSPOSE',m,n - m,-one,b(1,m + p - n + 1),ldb,y(m + p - n + 1),1, &
                     one,d,1)
           ! solve triangular system: r11*x = d1
           if (m > 0) then
              call la_strtrs('UPPER','NO TRANSPOSE','NON UNIT',m,1,a,lda,d,m,info)

              if (info > 0) then
                 info = 2
                 return
              end if
              ! copy d to x
              call la_scopy(m,d,1,x,1)
           end if
           ! backward transformation y = z**t *y
           call la_sormrq('LEFT','TRANSPOSE',p,1,np,b(max(1,n - p + 1),1),ldb,work( &
                     m + 1),y,max(1,p),work(m + np + 1),lwork - m - np,info)
           work(1) = m + np + max(lopt,int(work(m + np + 1),KIND=ilp))
           return
     end subroutine la_sggglm
     !> DGGGLM: solves a general Gauss-Markov linear model (GLM) problem:
     !> minimize || y ||_2   subject to   d = A*x + B*y
     !> x
     !> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     !> given N-vector. It is assumed that M <= N <= M+P, and
     !> rank(A) = M    and    rank( A B ) = N.
     !> Under these assumptions, the constrained equation is always
     !> consistent, and there is a unique solution x and a minimal 2-norm
     !> solution y, which is obtained using a generalized QR factorization
     !> of the matrices (A, B) given by
     !> A = Q*(R),   B = Q*T*Z.
     !> (0)
     !> In particular, if matrix B is square nonsingular, then the problem
     !> GLM is equivalent to the following weighted linear least squares
     !> problem
     !> minimize || inv(B)*(d-A*x) ||_2
     !> x
     !> where inv(B) denotes the inverse of B.

     pure subroutine la_dggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*),d(*)
           real(dp),intent(out) :: work(*),x(*),y(*)
        ! ===================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,lopt,lwkmin,lwkopt,nb,nb1,nb2,nb3,nb4,np
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           np = min(n,p)
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -2
           else if (p < 0 .or. p < n - m) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'DGEQRF',' ',n,m,-1,-1)
                 nb2 = la_ilaenv(1,'DGERQF',' ',n,m,-1,-1)
                 nb3 = la_ilaenv(1,'DORMQR',' ',n,m,p,-1)
                 nb4 = la_ilaenv(1,'DORMRQ',' ',n,m,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = m + np + max(n,p)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGGGLM',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              do i = 1,m
                 x(i) = zero
              end do
              do i = 1,p
                 y(i) = zero
              end do
              return
           end if
           ! compute the gqr factorization of matrices a and b:
                ! q**t*a = ( r11 ) m,    q**t*b*z**t = ( t11   t12 ) m
                         ! (  0  ) n-m                 (  0    t22 ) n-m
                            ! m                         m+p-n  n-m
           ! where r11 and t22 are upper triangular, and q and z are
           ! orthogonal.
           call la_dggqrf(n,m,p,a,lda,work,b,ldb,work(m + 1),work(m + np + 1),lwork - m - &
                     np,info)
           lopt = work(m + np + 1)
           ! update left-hand-side vector d = q**t*d = ( d1 ) m
                                                     ! ( d2 ) n-m
           call la_dormqr('LEFT','TRANSPOSE',n,1,m,a,lda,work,d,max(1,n),work(m + &
                     np + 1),lwork - m - np,info)
           lopt = max(lopt,int(work(m + np + 1),KIND=ilp))
           ! solve t22*y2 = d2 for y2
           if (n > m) then
              call la_dtrtrs('UPPER','NO TRANSPOSE','NON UNIT',n - m,1,b(m + 1,m + p - n + 1), &
                        ldb,d(m + 1),n - m,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              call la_dcopy(n - m,d(m + 1),1,y(m + p - n + 1),1)
           end if
           ! set y1 = 0
           do i = 1,m + p - n
              y(i) = zero
           end do
           ! update d1 = d1 - t12*y2
           call la_dgemv('NO TRANSPOSE',m,n - m,-one,b(1,m + p - n + 1),ldb,y(m + p - n + 1),1, &
                     one,d,1)
           ! solve triangular system: r11*x = d1
           if (m > 0) then
              call la_dtrtrs('UPPER','NO TRANSPOSE','NON UNIT',m,1,a,lda,d,m,info)

              if (info > 0) then
                 info = 2
                 return
              end if
              ! copy d to x
              call la_dcopy(m,d,1,x,1)
           end if
           ! backward transformation y = z**t *y
           call la_dormrq('LEFT','TRANSPOSE',p,1,np,b(max(1,n - p + 1),1),ldb,work( &
                     m + 1),y,max(1,p),work(m + np + 1),lwork - m - np,info)
           work(1) = m + np + max(lopt,int(work(m + np + 1),KIND=ilp))
           return
     end subroutine la_dggglm
     !> QGGGLM: solves a general Gauss-Markov linear model (GLM) problem:
     !> minimize || y ||_2   subject to   d = A*x + B*y
     !> x
     !> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     !> given N-vector. It is assumed that M <= N <= M+P, and
     !> rank(A) = M    and    rank( A B ) = N.
     !> Under these assumptions, the constrained equation is always
     !> consistent, and there is a unique solution x and a minimal 2-norm
     !> solution y, which is obtained using a generalized QR factorization
     !> of the matrices (A, B) given by
     !> A = Q*(R),   B = Q*T*Z.
     !> (0)
     !> In particular, if matrix B is square nonsingular, then the problem
     !> GLM is equivalent to the following weighted linear least squares
     !> problem
     !> minimize || inv(B)*(d-A*x) ||_2
     !> x
     !> where inv(B) denotes the inverse of B.

     pure subroutine la_qggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*),d(*)
           real(qp),intent(out) :: work(*),x(*),y(*)
        ! ===================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,lopt,lwkmin,lwkopt,nb,nb1,nb2,nb3,nb4,np
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           np = min(n,p)
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -2
           else if (p < 0 .or. p < n - m) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'QGEQRF',' ',n,m,-1,-1)
                 nb2 = la_ilaenv(1,'QGERQF',' ',n,m,-1,-1)
                 nb3 = la_ilaenv(1,'QORMQR',' ',n,m,p,-1)
                 nb4 = la_ilaenv(1,'QORMRQ',' ',n,m,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = m + np + max(n,p)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGGGLM',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              do i = 1,m
                 x(i) = zero
              end do
              do i = 1,p
                 y(i) = zero
              end do
              return
           end if
           ! compute the gqr factorization of matrices a and b:
                ! q**t*a = ( r11 ) m,    q**t*b*z**t = ( t11   t12 ) m
                         ! (  0  ) n-m                 (  0    t22 ) n-m
                            ! m                         m+p-n  n-m
           ! where r11 and t22 are upper triangular, and q and z are
           ! orthogonal.
           call la_qggqrf(n,m,p,a,lda,work,b,ldb,work(m + 1),work(m + np + 1),lwork - m - &
                     np,info)
           lopt = work(m + np + 1)
           ! update left-hand-side vector d = q**t*d = ( d1 ) m
                                                     ! ( d2 ) n-m
           call la_qormqr('LEFT','TRANSPOSE',n,1,m,a,lda,work,d,max(1,n),work(m + &
                     np + 1),lwork - m - np,info)
           lopt = max(lopt,int(work(m + np + 1),KIND=ilp))
           ! solve t22*y2 = d2 for y2
           if (n > m) then
              call la_qtrtrs('UPPER','NO TRANSPOSE','NON UNIT',n - m,1,b(m + 1,m + p - n + 1), &
                        ldb,d(m + 1),n - m,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              call la_qcopy(n - m,d(m + 1),1,y(m + p - n + 1),1)
           end if
           ! set y1 = 0
           do i = 1,m + p - n
              y(i) = zero
           end do
           ! update d1 = d1 - t12*y2
           call la_qgemv('NO TRANSPOSE',m,n - m,-one,b(1,m + p - n + 1),ldb,y(m + p - n + 1),1, &
                     one,d,1)
           ! solve triangular system: r11*x = d1
           if (m > 0) then
              call la_qtrtrs('UPPER','NO TRANSPOSE','NON UNIT',m,1,a,lda,d,m,info)

              if (info > 0) then
                 info = 2
                 return
              end if
              ! copy d to x
              call la_qcopy(m,d,1,x,1)
           end if
           ! backward transformation y = z**t *y
           call la_qormrq('LEFT','TRANSPOSE',p,1,np,b(max(1,n - p + 1),1),ldb,work( &
                     m + 1),y,max(1,p),work(m + np + 1),lwork - m - np,info)
           work(1) = m + np + max(lopt,int(work(m + np + 1),KIND=ilp))
           return
     end subroutine la_qggglm

     !> SGGLSE: solves the linear equality-constrained least squares (LSE)
     !> problem:
     !> minimize || c - A*x ||_2   subject to   B*x = d
     !> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     !> M-vector, and d is a given P-vector. It is assumed that
     !> P <= N <= M+P, and
     !> rank(B) = P and  rank( (A) ) = N.
     !> ( (B) )
     !> These conditions ensure that the LSE problem has a unique solution,
     !> which is obtained using a generalized RQ factorization of the
     !> matrices (B, A) given by
     !> B = (0 R)*Q,   A = Z*T*Q.

     pure subroutine la_sgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*),c(*),d(*)
           real(sp),intent(out) :: work(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lopt,lwkmin,lwkopt,mn,nb,nb1,nb2,nb3,nb4,nr
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (p < 0 .or. p > n .or. p < n - m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,p)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'SGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'SGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'SORMQR',' ',m,n,p,-1)
                 nb4 = la_ilaenv(1,'SORMRQ',' ',m,n,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = p + mn + max(m,n)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGGLSE',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! compute the grq factorization of matrices b and a:
                  ! b*q**t = (  0  t12 ) p   z**t*a*q**t = ( r11 r12 ) n-p
                              ! n-p  p                     (  0  r22 ) m+p-n
                                                            ! n-p  p
           ! where t12 and r11 are upper triangular, and q and z are
           ! orthogonal.
           call la_sggrqf(p,m,n,b,ldb,work,a,lda,work(p + 1),work(p + mn + 1),lwork - p - &
                     mn,info)
           lopt = work(p + mn + 1)
           ! update c = z**t *c = ( c1 ) n-p
                                ! ( c2 ) m+p-n
           call la_sormqr('LEFT','TRANSPOSE',m,1,mn,a,lda,work(p + 1),c,max(1,m), &
                     work(p + mn + 1),lwork - p - mn,info)
           lopt = max(lopt,int(work(p + mn + 1),KIND=ilp))
           ! solve t12*x2 = d for x2
           if (p > 0) then
              call la_strtrs('UPPER','NO TRANSPOSE','NON-UNIT',p,1,b(1,n - p + 1),ldb,d, &
                         p,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              ! put the solution in x
              call la_scopy(p,d,1,x(n - p + 1),1)
              ! update c1
              call la_sgemv('NO TRANSPOSE',n - p,p,-one,a(1,n - p + 1),lda,d,1,one,c,1 &
                        )
           end if
           ! solve r11*x1 = c1 for x1
           if (n > p) then
              call la_strtrs('UPPER','NO TRANSPOSE','NON-UNIT',n - p,1,a,lda,c,n - p, &
                        info)
              if (info > 0) then
                 info = 2
                 return
              end if
              ! put the solutions in x
              call la_scopy(n - p,c,1,x,1)
           end if
           ! compute the residual vector:
           if (m < n) then
              nr = m + p - n
              if (nr > 0) call la_sgemv('NO TRANSPOSE',nr,n - m,-one,a(n - p + 1,m + 1),lda,d( &
                        nr + 1),1,one,c(n - p + 1),1)
           else
              nr = p
           end if
           if (nr > 0) then
              call la_strmv('UPPER','NO TRANSPOSE','NON UNIT',nr,a(n - p + 1,n - p + 1),lda, &
                        d,1)
              call la_saxpy(nr,-one,d,1,c(n - p + 1),1)
           end if
           ! backward transformation x = q**t*x
           call la_sormrq('LEFT','TRANSPOSE',n,1,p,b,ldb,work(1),x,n,work(p + mn + 1 &
                     ),lwork - p - mn,info)
           work(1) = p + mn + max(lopt,int(work(p + mn + 1),KIND=ilp))
           return
     end subroutine la_sgglse
     !> DGGLSE: solves the linear equality-constrained least squares (LSE)
     !> problem:
     !> minimize || c - A*x ||_2   subject to   B*x = d
     !> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     !> M-vector, and d is a given P-vector. It is assumed that
     !> P <= N <= M+P, and
     !> rank(B) = P and  rank( (A) ) = N.
     !> ( (B) )
     !> These conditions ensure that the LSE problem has a unique solution,
     !> which is obtained using a generalized RQ factorization of the
     !> matrices (B, A) given by
     !> B = (0 R)*Q,   A = Z*T*Q.

     pure subroutine la_dgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*),c(*),d(*)
           real(dp),intent(out) :: work(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lopt,lwkmin,lwkopt,mn,nb,nb1,nb2,nb3,nb4,nr
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (p < 0 .or. p > n .or. p < n - m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,p)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'DGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'DGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'DORMQR',' ',m,n,p,-1)
                 nb4 = la_ilaenv(1,'DORMRQ',' ',m,n,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = p + mn + max(m,n)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGGLSE',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! compute the grq factorization of matrices b and a:
                  ! b*q**t = (  0  t12 ) p   z**t*a*q**t = ( r11 r12 ) n-p
                              ! n-p  p                     (  0  r22 ) m+p-n
                                                            ! n-p  p
           ! where t12 and r11 are upper triangular, and q and z are
           ! orthogonal.
           call la_dggrqf(p,m,n,b,ldb,work,a,lda,work(p + 1),work(p + mn + 1),lwork - p - &
                     mn,info)
           lopt = work(p + mn + 1)
           ! update c = z**t *c = ( c1 ) n-p
                                ! ( c2 ) m+p-n
           call la_dormqr('LEFT','TRANSPOSE',m,1,mn,a,lda,work(p + 1),c,max(1,m), &
                     work(p + mn + 1),lwork - p - mn,info)
           lopt = max(lopt,int(work(p + mn + 1),KIND=ilp))
           ! solve t12*x2 = d for x2
           if (p > 0) then
              call la_dtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',p,1,b(1,n - p + 1),ldb,d, &
                         p,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              ! put the solution in x
              call la_dcopy(p,d,1,x(n - p + 1),1)
              ! update c1
              call la_dgemv('NO TRANSPOSE',n - p,p,-one,a(1,n - p + 1),lda,d,1,one,c,1 &
                        )
           end if
           ! solve r11*x1 = c1 for x1
           if (n > p) then
              call la_dtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n - p,1,a,lda,c,n - p, &
                        info)
              if (info > 0) then
                 info = 2
                 return
              end if
              ! put the solutions in x
              call la_dcopy(n - p,c,1,x,1)
           end if
           ! compute the residual vector:
           if (m < n) then
              nr = m + p - n
              if (nr > 0) call la_dgemv('NO TRANSPOSE',nr,n - m,-one,a(n - p + 1,m + 1),lda,d( &
                        nr + 1),1,one,c(n - p + 1),1)
           else
              nr = p
           end if
           if (nr > 0) then
              call la_dtrmv('UPPER','NO TRANSPOSE','NON UNIT',nr,a(n - p + 1,n - p + 1),lda, &
                        d,1)
              call la_daxpy(nr,-one,d,1,c(n - p + 1),1)
           end if
           ! backward transformation x = q**t*x
           call la_dormrq('LEFT','TRANSPOSE',n,1,p,b,ldb,work(1),x,n,work(p + mn + 1 &
                     ),lwork - p - mn,info)
           work(1) = p + mn + max(lopt,int(work(p + mn + 1),KIND=ilp))
           return
     end subroutine la_dgglse
     !> QGGLSE: solves the linear equality-constrained least squares (LSE)
     !> problem:
     !> minimize || c - A*x ||_2   subject to   B*x = d
     !> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     !> M-vector, and d is a given P-vector. It is assumed that
     !> P <= N <= M+P, and
     !> rank(B) = P and  rank( (A) ) = N.
     !> ( (B) )
     !> These conditions ensure that the LSE problem has a unique solution,
     !> which is obtained using a generalized RQ factorization of the
     !> matrices (B, A) given by
     !> B = (0 R)*Q,   A = Z*T*Q.

     pure subroutine la_qgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*),c(*),d(*)
           real(qp),intent(out) :: work(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lopt,lwkmin,lwkopt,mn,nb,nb1,nb2,nb3,nb4,nr
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (p < 0 .or. p > n .or. p < n - m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,p)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'QGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'QGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'QORMQR',' ',m,n,p,-1)
                 nb4 = la_ilaenv(1,'QORMRQ',' ',m,n,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = p + mn + max(m,n)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGGLSE',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! compute the grq factorization of matrices b and a:
                  ! b*q**t = (  0  t12 ) p   z**t*a*q**t = ( r11 r12 ) n-p
                              ! n-p  p                     (  0  r22 ) m+p-n
                                                            ! n-p  p
           ! where t12 and r11 are upper triangular, and q and z are
           ! orthogonal.
           call la_qggrqf(p,m,n,b,ldb,work,a,lda,work(p + 1),work(p + mn + 1),lwork - p - &
                     mn,info)
           lopt = work(p + mn + 1)
           ! update c = z**t *c = ( c1 ) n-p
                                ! ( c2 ) m+p-n
           call la_qormqr('LEFT','TRANSPOSE',m,1,mn,a,lda,work(p + 1),c,max(1,m), &
                     work(p + mn + 1),lwork - p - mn,info)
           lopt = max(lopt,int(work(p + mn + 1),KIND=ilp))
           ! solve t12*x2 = d for x2
           if (p > 0) then
              call la_qtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',p,1,b(1,n - p + 1),ldb,d, &
                         p,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              ! put the solution in x
              call la_qcopy(p,d,1,x(n - p + 1),1)
              ! update c1
              call la_qgemv('NO TRANSPOSE',n - p,p,-one,a(1,n - p + 1),lda,d,1,one,c,1 &
                        )
           end if
           ! solve r11*x1 = c1 for x1
           if (n > p) then
              call la_qtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n - p,1,a,lda,c,n - p, &
                        info)
              if (info > 0) then
                 info = 2
                 return
              end if
              ! put the solutions in x
              call la_qcopy(n - p,c,1,x,1)
           end if
           ! compute the residual vector:
           if (m < n) then
              nr = m + p - n
              if (nr > 0) call la_qgemv('NO TRANSPOSE',nr,n - m,-one,a(n - p + 1,m + 1),lda,d( &
                        nr + 1),1,one,c(n - p + 1),1)
           else
              nr = p
           end if
           if (nr > 0) then
              call la_qtrmv('UPPER','NO TRANSPOSE','NON UNIT',nr,a(n - p + 1,n - p + 1),lda, &
                        d,1)
              call la_qaxpy(nr,-one,d,1,c(n - p + 1),1)
           end if
           ! backward transformation x = q**t*x
           call la_qormrq('LEFT','TRANSPOSE',n,1,p,b,ldb,work(1),x,n,work(p + mn + 1 &
                     ),lwork - p - mn,info)
           work(1) = p + mn + max(lopt,int(work(p + mn + 1),KIND=ilp))
           return
     end subroutine la_qgglse

     !> CGGGLM: solves a general Gauss-Markov linear model (GLM) problem:
     !> minimize || y ||_2   subject to   d = A*x + B*y
     !> x
     !> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     !> given N-vector. It is assumed that M <= N <= M+P, and
     !> rank(A) = M    and    rank( A B ) = N.
     !> Under these assumptions, the constrained equation is always
     !> consistent, and there is a unique solution x and a minimal 2-norm
     !> solution y, which is obtained using a generalized QR factorization
     !> of the matrices (A, B) given by
     !> A = Q*(R),   B = Q*T*Z.
     !> (0)
     !> In particular, if matrix B is square nonsingular, then the problem
     !> GLM is equivalent to the following weighted linear least squares
     !> problem
     !> minimize || inv(B)*(d-A*x) ||_2
     !> x
     !> where inv(B) denotes the inverse of B.

     pure subroutine la_cggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*),d(*)
           complex(sp),intent(out) :: work(*),x(*),y(*)
        ! ===================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,lopt,lwkmin,lwkopt,nb,nb1,nb2,nb3,nb4,np
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           np = min(n,p)
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -2
           else if (p < 0 .or. p < n - m) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'CGEQRF',' ',n,m,-1,-1)
                 nb2 = la_ilaenv(1,'CGERQF',' ',n,m,-1,-1)
                 nb3 = la_ilaenv(1,'CUNMQR',' ',n,m,p,-1)
                 nb4 = la_ilaenv(1,'CUNMRQ',' ',n,m,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = m + np + max(n,p)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGGGLM',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              do i = 1,m
                 x(i) = czero
              end do
              do i = 1,p
                 y(i) = czero
              end do
              return
           end if
           ! compute the gqr factorization of matrices a and b:
                ! q**h*a = ( r11 ) m,    q**h*b*z**h = ( t11   t12 ) m
                         ! (  0  ) n-m                 (  0    t22 ) n-m
                            ! m                         m+p-n  n-m
           ! where r11 and t22 are upper triangular, and q and z are
           ! unitary.
           call la_cggqrf(n,m,p,a,lda,work,b,ldb,work(m + 1),work(m + np + 1),lwork - m - &
                     np,info)
           lopt = real(work(m + np + 1),KIND=sp)
           ! update left-hand-side vector d = q**h*d = ( d1 ) m
                                                     ! ( d2 ) n-m
           call la_cunmqr('LEFT','CONJUGATE TRANSPOSE',n,1,m,a,lda,work,d,max(1,n) &
                     ,work(m + np + 1),lwork - m - np,info)
           lopt = max(lopt,int(work(m + np + 1),KIND=ilp))
           ! solve t22*y2 = d2 for y2
           if (n > m) then
              call la_ctrtrs('UPPER','NO TRANSPOSE','NON UNIT',n - m,1,b(m + 1,m + p - n + 1), &
                        ldb,d(m + 1),n - m,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              call la_ccopy(n - m,d(m + 1),1,y(m + p - n + 1),1)
           end if
           ! set y1 = 0
           do i = 1,m + p - n
              y(i) = czero
           end do
           ! update d1 = d1 - t12*y2
           call la_cgemv('NO TRANSPOSE',m,n - m,-cone,b(1,m + p - n + 1),ldb,y(m + p - n + 1),1, &
                      cone,d,1)
           ! solve triangular system: r11*x = d1
           if (m > 0) then
              call la_ctrtrs('UPPER','NO TRANSPOSE','NON UNIT',m,1,a,lda,d,m,info)

              if (info > 0) then
                 info = 2
                 return
              end if
              ! copy d to x
              call la_ccopy(m,d,1,x,1)
           end if
           ! backward transformation y = z**h *y
           call la_cunmrq('LEFT','CONJUGATE TRANSPOSE',p,1,np,b(max(1,n - p + 1),1), &
                     ldb,work(m + 1),y,max(1,p),work(m + np + 1),lwork - m - np,info)
           work(1) = m + np + max(lopt,int(work(m + np + 1),KIND=ilp))
           return
     end subroutine la_cggglm
     !> ZGGGLM: solves a general Gauss-Markov linear model (GLM) problem:
     !> minimize || y ||_2   subject to   d = A*x + B*y
     !> x
     !> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     !> given N-vector. It is assumed that M <= N <= M+P, and
     !> rank(A) = M    and    rank( A B ) = N.
     !> Under these assumptions, the constrained equation is always
     !> consistent, and there is a unique solution x and a minimal 2-norm
     !> solution y, which is obtained using a generalized QR factorization
     !> of the matrices (A, B) given by
     !> A = Q*(R),   B = Q*T*Z.
     !> (0)
     !> In particular, if matrix B is square nonsingular, then the problem
     !> GLM is equivalent to the following weighted linear least squares
     !> problem
     !> minimize || inv(B)*(d-A*x) ||_2
     !> x
     !> where inv(B) denotes the inverse of B.

     pure subroutine la_zggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*),d(*)
           complex(dp),intent(out) :: work(*),x(*),y(*)
        ! ===================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,lopt,lwkmin,lwkopt,nb,nb1,nb2,nb3,nb4,np
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           np = min(n,p)
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -2
           else if (p < 0 .or. p < n - m) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'ZGEQRF',' ',n,m,-1,-1)
                 nb2 = la_ilaenv(1,'ZGERQF',' ',n,m,-1,-1)
                 nb3 = la_ilaenv(1,'ZUNMQR',' ',n,m,p,-1)
                 nb4 = la_ilaenv(1,'ZUNMRQ',' ',n,m,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = m + np + max(n,p)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGGGLM',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              do i = 1,m
                 x(i) = czero
              end do
              do i = 1,p
                 y(i) = czero
              end do
              return
           end if
           ! compute the gqr factorization of matrices a and b:
                ! q**h*a = ( r11 ) m,    q**h*b*z**h = ( t11   t12 ) m
                         ! (  0  ) n-m                 (  0    t22 ) n-m
                            ! m                         m+p-n  n-m
           ! where r11 and t22 are upper triangular, and q and z are
           ! unitary.
           call la_zggqrf(n,m,p,a,lda,work,b,ldb,work(m + 1),work(m + np + 1),lwork - m - &
                     np,info)
           lopt = real(work(m + np + 1),KIND=dp)
           ! update left-hand-side vector d = q**h*d = ( d1 ) m
                                                     ! ( d2 ) n-m
           call la_zunmqr('LEFT','CONJUGATE TRANSPOSE',n,1,m,a,lda,work,d,max(1,n) &
                     ,work(m + np + 1),lwork - m - np,info)
           lopt = max(lopt,int(work(m + np + 1),KIND=ilp))
           ! solve t22*y2 = d2 for y2
           if (n > m) then
              call la_ztrtrs('UPPER','NO TRANSPOSE','NON UNIT',n - m,1,b(m + 1,m + p - n + 1), &
                        ldb,d(m + 1),n - m,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              call la_zcopy(n - m,d(m + 1),1,y(m + p - n + 1),1)
           end if
           ! set y1 = 0
           do i = 1,m + p - n
              y(i) = czero
           end do
           ! update d1 = d1 - t12*y2
           call la_zgemv('NO TRANSPOSE',m,n - m,-cone,b(1,m + p - n + 1),ldb,y(m + p - n + 1),1, &
                      cone,d,1)
           ! solve triangular system: r11*x = d1
           if (m > 0) then
              call la_ztrtrs('UPPER','NO TRANSPOSE','NON UNIT',m,1,a,lda,d,m,info)

              if (info > 0) then
                 info = 2
                 return
              end if
              ! copy d to x
              call la_zcopy(m,d,1,x,1)
           end if
           ! backward transformation y = z**h *y
           call la_zunmrq('LEFT','CONJUGATE TRANSPOSE',p,1,np,b(max(1,n - p + 1),1), &
                     ldb,work(m + 1),y,max(1,p),work(m + np + 1),lwork - m - np,info)
           work(1) = m + np + max(lopt,int(work(m + np + 1),KIND=ilp))
           return
     end subroutine la_zggglm
     !> WGGGLM: solves a general Gauss-Markov linear model (GLM) problem:
     !> minimize || y ||_2   subject to   d = A*x + B*y
     !> x
     !> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     !> given N-vector. It is assumed that M <= N <= M+P, and
     !> rank(A) = M    and    rank( A B ) = N.
     !> Under these assumptions, the constrained equation is always
     !> consistent, and there is a unique solution x and a minimal 2-norm
     !> solution y, which is obtained using a generalized QR factorization
     !> of the matrices (A, B) given by
     !> A = Q*(R),   B = Q*T*Z.
     !> (0)
     !> In particular, if matrix B is square nonsingular, then the problem
     !> GLM is equivalent to the following weighted linear least squares
     !> problem
     !> minimize || inv(B)*(d-A*x) ||_2
     !> x
     !> where inv(B) denotes the inverse of B.

     pure subroutine la_wggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*),d(*)
           complex(qp),intent(out) :: work(*),x(*),y(*)
        ! ===================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,lopt,lwkmin,lwkopt,nb,nb1,nb2,nb3,nb4,np
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           np = min(n,p)
           lquery = (lwork == -1)
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -2
           else if (p < 0 .or. p < n - m) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'WGEQRF',' ',n,m,-1,-1)
                 nb2 = la_ilaenv(1,'WGERQF',' ',n,m,-1,-1)
                 nb3 = la_ilaenv(1,'WUNMQR',' ',n,m,p,-1)
                 nb4 = la_ilaenv(1,'WUNMRQ',' ',n,m,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = m + np + max(n,p)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGGGLM',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              do i = 1,m
                 x(i) = czero
              end do
              do i = 1,p
                 y(i) = czero
              end do
              return
           end if
           ! compute the gqr factorization of matrices a and b:
                ! q**h*a = ( r11 ) m,    q**h*b*z**h = ( t11   t12 ) m
                         ! (  0  ) n-m                 (  0    t22 ) n-m
                            ! m                         m+p-n  n-m
           ! where r11 and t22 are upper triangular, and q and z are
           ! unitary.
           call la_wggqrf(n,m,p,a,lda,work,b,ldb,work(m + 1),work(m + np + 1),lwork - m - &
                     np,info)
           lopt = real(work(m + np + 1),KIND=qp)
           ! update left-hand-side vector d = q**h*d = ( d1 ) m
                                                     ! ( d2 ) n-m
           call la_wunmqr('LEFT','CONJUGATE TRANSPOSE',n,1,m,a,lda,work,d,max(1,n) &
                     ,work(m + np + 1),lwork - m - np,info)
           lopt = max(lopt,int(work(m + np + 1),KIND=ilp))
           ! solve t22*y2 = d2 for y2
           if (n > m) then
              call la_wtrtrs('UPPER','NO TRANSPOSE','NON UNIT',n - m,1,b(m + 1,m + p - n + 1), &
                        ldb,d(m + 1),n - m,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              call la_wcopy(n - m,d(m + 1),1,y(m + p - n + 1),1)
           end if
           ! set y1 = 0
           do i = 1,m + p - n
              y(i) = czero
           end do
           ! update d1 = d1 - t12*y2
           call la_wgemv('NO TRANSPOSE',m,n - m,-cone,b(1,m + p - n + 1),ldb,y(m + p - n + 1),1, &
                      cone,d,1)
           ! solve triangular system: r11*x = d1
           if (m > 0) then
              call la_wtrtrs('UPPER','NO TRANSPOSE','NON UNIT',m,1,a,lda,d,m,info)

              if (info > 0) then
                 info = 2
                 return
              end if
              ! copy d to x
              call la_wcopy(m,d,1,x,1)
           end if
           ! backward transformation y = z**h *y
           call la_wunmrq('LEFT','CONJUGATE TRANSPOSE',p,1,np,b(max(1,n - p + 1),1), &
                     ldb,work(m + 1),y,max(1,p),work(m + np + 1),lwork - m - np,info)
           work(1) = m + np + max(lopt,int(work(m + np + 1),KIND=ilp))
           return
     end subroutine la_wggglm

     !> CGGLSE: solves the linear equality-constrained least squares (LSE)
     !> problem:
     !> minimize || c - A*x ||_2   subject to   B*x = d
     !> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     !> M-vector, and d is a given P-vector. It is assumed that
     !> P <= N <= M+P, and
     !> rank(B) = P and  rank( (A) ) = N.
     !> ( (B) )
     !> These conditions ensure that the LSE problem has a unique solution,
     !> which is obtained using a generalized RQ factorization of the
     !> matrices (B, A) given by
     !> B = (0 R)*Q,   A = Z*T*Q.

     pure subroutine la_cgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*),c(*),d(*)
           complex(sp),intent(out) :: work(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lopt,lwkmin,lwkopt,mn,nb,nb1,nb2,nb3,nb4,nr
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (p < 0 .or. p > n .or. p < n - m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,p)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'CGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'CGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'CUNMQR',' ',m,n,p,-1)
                 nb4 = la_ilaenv(1,'CUNMRQ',' ',m,n,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = p + mn + max(m,n)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGGLSE',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! compute the grq factorization of matrices b and a:
                  ! b*q**h = (  0  t12 ) p   z**h*a*q**h = ( r11 r12 ) n-p
                              ! n-p  p                     (  0  r22 ) m+p-n
                                                            ! n-p  p
           ! where t12 and r11 are upper triangular, and q and z are
           ! unitary.
           call la_cggrqf(p,m,n,b,ldb,work,a,lda,work(p + 1),work(p + mn + 1),lwork - p - &
                     mn,info)
           lopt = real(work(p + mn + 1),KIND=sp)
           ! update c = z**h *c = ( c1 ) n-p
                             ! ( c2 ) m+p-n
           call la_cunmqr('LEFT','CONJUGATE TRANSPOSE',m,1,mn,a,lda,work(p + 1),c, &
                     max(1,m),work(p + mn + 1),lwork - p - mn,info)
           lopt = max(lopt,int(work(p + mn + 1),KIND=ilp))
           ! solve t12*x2 = d for x2
           if (p > 0) then
              call la_ctrtrs('UPPER','NO TRANSPOSE','NON-UNIT',p,1,b(1,n - p + 1),ldb,d, &
                         p,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              ! put the solution in x
              call la_ccopy(p,d,1,x(n - p + 1),1)
              ! update c1
              call la_cgemv('NO TRANSPOSE',n - p,p,-cone,a(1,n - p + 1),lda,d,1,cone,c, &
                        1)
           end if
           ! solve r11*x1 = c1 for x1
           if (n > p) then
              call la_ctrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n - p,1,a,lda,c,n - p, &
                        info)
              if (info > 0) then
                 info = 2
                 return
              end if
              ! put the solutions in x
              call la_ccopy(n - p,c,1,x,1)
           end if
           ! compute the residual vector:
           if (m < n) then
              nr = m + p - n
              if (nr > 0) call la_cgemv('NO TRANSPOSE',nr,n - m,-cone,a(n - p + 1,m + 1),lda,d( &
                         nr + 1),1,cone,c(n - p + 1),1)
           else
              nr = p
           end if
           if (nr > 0) then
              call la_ctrmv('UPPER','NO TRANSPOSE','NON UNIT',nr,a(n - p + 1,n - p + 1),lda, &
                        d,1)
              call la_caxpy(nr,-cone,d,1,c(n - p + 1),1)
           end if
           ! backward transformation x = q**h*x
           call la_cunmrq('LEFT','CONJUGATE TRANSPOSE',n,1,p,b,ldb,work(1),x,n, &
                     work(p + mn + 1),lwork - p - mn,info)
           work(1) = p + mn + max(lopt,int(work(p + mn + 1),KIND=ilp))
           return
     end subroutine la_cgglse
     !> ZGGLSE: solves the linear equality-constrained least squares (LSE)
     !> problem:
     !> minimize || c - A*x ||_2   subject to   B*x = d
     !> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     !> M-vector, and d is a given P-vector. It is assumed that
     !> P <= N <= M+P, and
     !> rank(B) = P and  rank( (A) ) = N.
     !> ( (B) )
     !> These conditions ensure that the LSE problem has a unique solution,
     !> which is obtained using a generalized RQ factorization of the
     !> matrices (B, A) given by
     !> B = (0 R)*Q,   A = Z*T*Q.

     pure subroutine la_zgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*),c(*),d(*)
           complex(dp),intent(out) :: work(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lopt,lwkmin,lwkopt,mn,nb,nb1,nb2,nb3,nb4,nr
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (p < 0 .or. p > n .or. p < n - m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,p)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'ZGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'ZGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'ZUNMQR',' ',m,n,p,-1)
                 nb4 = la_ilaenv(1,'ZUNMRQ',' ',m,n,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = p + mn + max(m,n)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGGLSE',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! compute the grq factorization of matrices b and a:
                  ! b*q**h = (  0  t12 ) p   z**h*a*q**h = ( r11 r12 ) n-p
                              ! n-p  p                     (  0  r22 ) m+p-n
                                                            ! n-p  p
           ! where t12 and r11 are upper triangular, and q and z are
           ! unitary.
           call la_zggrqf(p,m,n,b,ldb,work,a,lda,work(p + 1),work(p + mn + 1),lwork - p - &
                     mn,info)
           lopt = real(work(p + mn + 1),KIND=dp)
           ! update c = z**h *c = ( c1 ) n-p
                             ! ( c2 ) m+p-n
           call la_zunmqr('LEFT','CONJUGATE TRANSPOSE',m,1,mn,a,lda,work(p + 1),c, &
                     max(1,m),work(p + mn + 1),lwork - p - mn,info)
           lopt = max(lopt,int(work(p + mn + 1),KIND=ilp))
           ! solve t12*x2 = d for x2
           if (p > 0) then
              call la_ztrtrs('UPPER','NO TRANSPOSE','NON-UNIT',p,1,b(1,n - p + 1),ldb,d, &
                         p,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              ! put the solution in x
              call la_zcopy(p,d,1,x(n - p + 1),1)
              ! update c1
              call la_zgemv('NO TRANSPOSE',n - p,p,-cone,a(1,n - p + 1),lda,d,1,cone,c, &
                        1)
           end if
           ! solve r11*x1 = c1 for x1
           if (n > p) then
              call la_ztrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n - p,1,a,lda,c,n - p, &
                        info)
              if (info > 0) then
                 info = 2
                 return
              end if
              ! put the solutions in x
              call la_zcopy(n - p,c,1,x,1)
           end if
           ! compute the residual vector:
           if (m < n) then
              nr = m + p - n
              if (nr > 0) call la_zgemv('NO TRANSPOSE',nr,n - m,-cone,a(n - p + 1,m + 1),lda,d( &
                         nr + 1),1,cone,c(n - p + 1),1)
           else
              nr = p
           end if
           if (nr > 0) then
              call la_ztrmv('UPPER','NO TRANSPOSE','NON UNIT',nr,a(n - p + 1,n - p + 1),lda, &
                        d,1)
              call la_zaxpy(nr,-cone,d,1,c(n - p + 1),1)
           end if
           ! backward transformation x = q**h*x
           call la_zunmrq('LEFT','CONJUGATE TRANSPOSE',n,1,p,b,ldb,work(1),x,n, &
                     work(p + mn + 1),lwork - p - mn,info)
           work(1) = p + mn + max(lopt,int(work(p + mn + 1),KIND=ilp))
           return
     end subroutine la_zgglse
     !> WGGLSE: solves the linear equality-constrained least squares (LSE)
     !> problem:
     !> minimize || c - A*x ||_2   subject to   B*x = d
     !> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     !> M-vector, and d is a given P-vector. It is assumed that
     !> P <= N <= M+P, and
     !> rank(B) = P and  rank( (A) ) = N.
     !> ( (B) )
     !> These conditions ensure that the LSE problem has a unique solution,
     !> which is obtained using a generalized RQ factorization of the
     !> matrices (B, A) given by
     !> B = (0 R)*Q,   A = Z*T*Q.

     pure subroutine la_wgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,m,n,p
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*),c(*),d(*)
           complex(qp),intent(out) :: work(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lopt,lwkmin,lwkopt,mn,nb,nb1,nb2,nb3,nb4,nr
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           mn = min(m,n)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (p < 0 .or. p > n .or. p < n - m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,p)) then
              info = -7
           end if
           ! calculate workspace
           if (info == 0) then
              if (n == 0) then
                 lwkmin = 1
                 lwkopt = 1
              else
                 nb1 = la_ilaenv(1,'WGEQRF',' ',m,n,-1,-1)
                 nb2 = la_ilaenv(1,'WGERQF',' ',m,n,-1,-1)
                 nb3 = la_ilaenv(1,'WUNMQR',' ',m,n,p,-1)
                 nb4 = la_ilaenv(1,'WUNMRQ',' ',m,n,p,-1)
                 nb = max(nb1,nb2,nb3,nb4)
                 lwkmin = m + n + p
                 lwkopt = p + mn + max(m,n)*nb
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGGLSE',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! compute the grq factorization of matrices b and a:
                  ! b*q**h = (  0  t12 ) p   z**h*a*q**h = ( r11 r12 ) n-p
                              ! n-p  p                     (  0  r22 ) m+p-n
                                                            ! n-p  p
           ! where t12 and r11 are upper triangular, and q and z are
           ! unitary.
           call la_wggrqf(p,m,n,b,ldb,work,a,lda,work(p + 1),work(p + mn + 1),lwork - p - &
                     mn,info)
           lopt = real(work(p + mn + 1),KIND=qp)
           ! update c = z**h *c = ( c1 ) n-p
                             ! ( c2 ) m+p-n
           call la_wunmqr('LEFT','CONJUGATE TRANSPOSE',m,1,mn,a,lda,work(p + 1),c, &
                     max(1,m),work(p + mn + 1),lwork - p - mn,info)
           lopt = max(lopt,int(work(p + mn + 1),KIND=ilp))
           ! solve t12*x2 = d for x2
           if (p > 0) then
              call la_wtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',p,1,b(1,n - p + 1),ldb,d, &
                         p,info)
              if (info > 0) then
                 info = 1
                 return
              end if
              ! put the solution in x
              call la_wcopy(p,d,1,x(n - p + 1),1)
              ! update c1
              call la_wgemv('NO TRANSPOSE',n - p,p,-cone,a(1,n - p + 1),lda,d,1,cone,c, &
                        1)
           end if
           ! solve r11*x1 = c1 for x1
           if (n > p) then
              call la_wtrtrs('UPPER','NO TRANSPOSE','NON-UNIT',n - p,1,a,lda,c,n - p, &
                        info)
              if (info > 0) then
                 info = 2
                 return
              end if
              ! put the solutions in x
              call la_wcopy(n - p,c,1,x,1)
           end if
           ! compute the residual vector:
           if (m < n) then
              nr = m + p - n
              if (nr > 0) call la_wgemv('NO TRANSPOSE',nr,n - m,-cone,a(n - p + 1,m + 1),lda,d( &
                         nr + 1),1,cone,c(n - p + 1),1)
           else
              nr = p
           end if
           if (nr > 0) then
              call la_wtrmv('UPPER','NO TRANSPOSE','NON UNIT',nr,a(n - p + 1,n - p + 1),lda, &
                        d,1)
              call la_waxpy(nr,-cone,d,1,c(n - p + 1),1)
           end if
           ! backward transformation x = q**h*x
           call la_wunmrq('LEFT','CONJUGATE TRANSPOSE',n,1,p,b,ldb,work(1),x,n, &
                     work(p + mn + 1),lwork - p - mn,info)
           work(1) = p + mn + max(lopt,int(work(p + mn + 1),KIND=ilp))
           return
     end subroutine la_wgglse

end module la_lapack_lsq_constrained
