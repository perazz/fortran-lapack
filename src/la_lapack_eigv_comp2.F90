!> Generalized nonsymmetric eigenproblem components: eigenvectors, block swaps, deflating subspaces, Sylvester solves
module la_lapack_eigv_comp2
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level3_gen
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_scalar
     use la_lapack_eigv_gen_aux
     use la_lapack_eigv_sym_comp
     use la_lapack_givens_jacobi_rot
     use la_lapack_orthogonal_factors_qr
     use la_lapack_solve_aux
     use la_lapack_solve_lu_comp
     use la_lapack_svd_comp2
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_stgsy2
     public :: la_stgsyl
     public :: la_slagv2
     public :: la_stgevc
     public :: la_stgex2
     public :: la_stgexc
     public :: la_stgsen
     public :: la_stgsna
     public :: la_dtgsy2
     public :: la_dtgsyl
     public :: la_dlagv2
     public :: la_dtgevc
     public :: la_dtgex2
     public :: la_dtgexc
     public :: la_dtgsen
     public :: la_dtgsna
#ifdef LA_WITH_XDP
     public :: la_xtgsy2
     public :: la_xtgsyl
     public :: la_xlagv2
     public :: la_xtgevc
     public :: la_xtgex2
     public :: la_xtgexc
     public :: la_xtgsen
     public :: la_xtgsna
#endif
#ifdef LA_WITH_QP
     public :: la_qtgsy2
     public :: la_qtgsyl
     public :: la_qlagv2
     public :: la_qtgevc
     public :: la_qtgex2
     public :: la_qtgexc
     public :: la_qtgsen
     public :: la_qtgsna
#endif
     public :: la_ctgevc
     public :: la_ctgex2
     public :: la_ctgexc
     public :: la_ctgsy2
     public :: la_ctgsyl
     public :: la_ctgsen
     public :: la_ctgsna
     public :: la_ztgevc
     public :: la_ztgex2
     public :: la_ztgexc
     public :: la_ztgsy2
     public :: la_ztgsyl
     public :: la_ztgsen
     public :: la_ztgsna
#ifdef LA_WITH_XDP
     public :: la_ytgevc
     public :: la_ytgex2
     public :: la_ytgexc
     public :: la_ytgsy2
     public :: la_ytgsyl
     public :: la_ytgsen
     public :: la_ytgsna
#endif
#ifdef LA_WITH_QP
     public :: la_wtgevc
     public :: la_wtgex2
     public :: la_wtgexc
     public :: la_wtgsy2
     public :: la_wtgsyl
     public :: la_wtgsen
     public :: la_wtgsna
#endif

     contains

     !> STGSY2: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                (1)
     !> D * R - L * E = scale * F,
     !> using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)
     !> must be in generalized Schur canonical form, i.e. A, B are upper
     !> quasi triangular and D, E are upper triangular. The solution (R, L)
     !> overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor
     !> chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Z*x = scale*b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**T, Im) ],
     !> Ik is the identity matrix of size k and X**T is the transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> In the process of solving (1), we solve a number of such systems
     !> where Dim(In), Dim(In) = 1 or 2.
     !> If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,
     !> which is equivalent to solve for R and L in
     !> A**T * R  + D**T * L   = scale * C           (3)
     !> R  * B**T + L  * E**T  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> sigma_min(Z) using reverse communication with SLACON.
     !> STGSY2 also (IJOB >= 1) contributes to the computation in STGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of the matrix pair in
     !> STGSYL. See STGSYL for details.

     pure subroutine la_stgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,iwork,pq,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info,pq
           real(sp),intent(inout) :: rdscal,rdsum
           real(sp),intent(out) :: scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(sp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
        ! replaced various illegal calls to la_scopy by calls to la_slaset.
        ! sven hammarling, 27/5/02.
           ! Parameters
           integer(ilp),parameter :: ldz = 8

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ie,ierr,ii,is,isp1,j,je,jj,js,jsp1,k,mb,nb,p,q, &
                     zdim
           real(sp) :: alpha,scaloc
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           real(sp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STGSY2',-info)
              return
           end if
           ! determine block structure of a
           pq = 0
           p = 0
           i = 1
           10 continue
           if (i > m) go to 20
           p = p + 1
           iwork(p) = i
           if (i == m) go to 20
           if (a(i + 1,i) /= zero) then
              i = i + 2
           else
              i = i + 1
           end if
           go to 10
           20 continue
           iwork(p + 1) = m + 1
           ! determine block structure of b
           q = p + 1
           j = 1
           30 continue
           if (j > n) go to 40
           q = q + 1
           iwork(q) = j
           if (j == n) go to 40
           if (b(j + 1,j) /= zero) then
              j = j + 2
           else
              j = j + 1
           end if
           go to 30
           40 continue
           iwork(q + 1) = n + 1
           pq = p*(q - p - 1)
           if (notran) then
              ! solve (i, j) - subsystem
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
              scale = one
              scaloc = one
              loop_120: do j = p + 2,q
                 js = iwork(j)
                 jsp1 = js + 1
                 je = iwork(j + 1) - 1
                 nb = je - js + 1
                 loop_110: do i = p,1,-1
                    is = iwork(i)
                    isp1 = is + 1
                    ie = iwork(i + 1) - 1
                    mb = ie - is + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = d(is,is)
                       z(1,2) = -b(js,js)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_sscal(m,scaloc,c(1,k),1)
                                call la_sscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_slatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          alpha = -rhs(1)
                          call la_saxpy(is - 1,alpha,a(1,is),1,c(1,js),1)
                          call la_saxpy(is - 1,alpha,d(1,is),1,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_saxpy(n - je,rhs(2),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_saxpy(n - je,rhs(2),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = zero
                       z(4,2) = d(is,is)
                       z(1,3) = -b(js,js)
                       z(2,3) = -b(js,jsp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = -e(js,jsp1)
                       z(1,4) = -b(jsp1,js)
                       z(2,4) = -b(jsp1,jsp1)
                       z(3,4) = zero
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_sscal(m,scaloc,c(1,k),1)
                                call la_sscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_slatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_sger(is - 1,nb,-one,a(1,is),1,rhs(1),1,c(1,js), &
                                     ldc)
                          call la_sger(is - 1,nb,-one,d(1,is),1,rhs(1),1,f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_saxpy(n - je,rhs(3),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_saxpy(n - je,rhs(3),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                          call la_saxpy(n - je,rhs(4),b(jsp1,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_saxpy(n - je,rhs(4),e(jsp1,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = d(is,isp1)
                       z(4,2) = d(isp1,isp1)
                       z(1,3) = -b(js,js)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = -b(js,js)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_sscal(m,scaloc,c(1,k),1)
                                call la_sscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_slatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_sgemv('N',is - 1,mb,-one,a(1,is),lda,rhs(1),1, &
                                    one,c(1,js),1)
                          call la_sgemv('N',is - 1,mb,-one,d(1,is),ldd,rhs(1),1, &
                                    one,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_sger(mb,n - je,one,rhs(3),1,b(js,je + 1),ldb,c(is, &
                                    je + 1),ldc)
                          call la_sger(mb,n - je,one,rhs(3),1,e(js,je + 1),lde,f(is, &
                                    je + 1),ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z * x = rhs
                       call la_slaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(5,1) = d(is,is)
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(5,2) = d(is,isp1)
                       z(6,2) = d(isp1,isp1)
                       z(3,3) = a(is,is)
                       z(4,3) = a(isp1,is)
                       z(7,3) = d(is,is)
                       z(3,4) = a(is,isp1)
                       z(4,4) = a(isp1,isp1)
                       z(7,4) = d(is,isp1)
                       z(8,4) = d(isp1,isp1)
                       z(1,5) = -b(js,js)
                       z(3,5) = -b(js,jsp1)
                       z(5,5) = -e(js,js)
                       z(7,5) = -e(js,jsp1)
                       z(2,6) = -b(js,js)
                       z(4,6) = -b(js,jsp1)
                       z(6,6) = -e(js,js)
                       z(8,6) = -e(js,jsp1)
                       z(1,7) = -b(jsp1,js)
                       z(3,7) = -b(jsp1,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(2,8) = -b(jsp1,js)
                       z(4,8) = -b(jsp1,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_scopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_scopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_sscal(m,scaloc,c(1,k),1)
                                call la_sscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_slatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_scopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_scopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_sgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,rhs(1 &
                                    ),mb,one,c(1,js),ldc)
                          call la_sgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,rhs(1 &
                                    ),mb,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          k = mb*nb + 1
                          call la_sgemm('N','N',mb,n - je,nb,one,rhs(k),mb,b(js,je + &
                                    1),ldb,one,c(is,je + 1),ldc)
                          call la_sgemm('N','N',mb,n - je,nb,one,rhs(k),mb,e(js,je + &
                                    1),lde,one,f(is,je + 1),ldf)
                       end if
                    end if
                 end do loop_110
              end do loop_120
           else
              ! solve (i, j) - subsystem
                   ! a(i, i)**t * r(i, j) + d(i, i)**t * l(j, j)  =  c(i, j)
                   ! r(i, i)  * b(j, j) + l(i, j)  * e(j, j)  = -f(i, j)
              ! for i = 1, 2, ..., p, j = q, q - 1, ..., 1
              scale = one
              scaloc = one
              loop_200: do i = 1,p
                 is = iwork(i)
                 isp1 = is + 1
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_190: do j = q,p + 2,-1
                    js = iwork(j)
                    jsp1 = js + 1
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = -b(js,js)
                       z(1,2) = d(is,is)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z**t * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_sscal(m,scaloc,c(1,k),1)
                             call la_sscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          alpha = rhs(1)
                          call la_saxpy(js - 1,alpha,b(1,js),1,f(is,1),ldf)
                          alpha = rhs(2)
                          call la_saxpy(js - 1,alpha,e(1,js),1,f(is,1),ldf)
                       end if
                       if (i < p) then
                          alpha = -rhs(1)
                          call la_saxpy(m - ie,alpha,a(is,ie + 1),lda,c(ie + 1,js),1)

                          alpha = -rhs(2)
                          call la_saxpy(m - ie,alpha,d(is,ie + 1),ldd,c(ie + 1,js),1)

                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = -b(js,js)
                       z(4,1) = -b(jsp1,js)
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = -b(js,jsp1)
                       z(4,2) = -b(jsp1,jsp1)
                       z(1,3) = d(is,is)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(is,is)
                       z(3,4) = -e(js,jsp1)
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z**t * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_sscal(m,scaloc,c(1,k),1)
                             call la_sscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_saxpy(js - 1,rhs(1),b(1,js),1,f(is,1),ldf)

                          call la_saxpy(js - 1,rhs(2),b(1,jsp1),1,f(is,1),ldf)

                          call la_saxpy(js - 1,rhs(3),e(1,js),1,f(is,1),ldf)

                          call la_saxpy(js - 1,rhs(4),e(1,jsp1),1,f(is,1),ldf)

                       end if
                       if (i < p) then
                          call la_sger(m - ie,nb,-one,a(is,ie + 1),lda,rhs(1),1,c(ie + &
                                    1,js),ldc)
                          call la_sger(m - ie,nb,-one,d(is,ie + 1),ldd,rhs(3),1,c(ie + &
                                    1,js),ldc)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(3,1) = -b(js,js)
                       z(4,1) = zero
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = zero
                       z(4,2) = -b(js,js)
                       z(1,3) = d(is,is)
                       z(2,3) = d(is,isp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(isp1,isp1)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z**t * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_sscal(m,scaloc,c(1,k),1)
                             call la_sscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_sger(mb,js - 1,one,rhs(1),1,b(1,js),1,f(is,1), &
                                    ldf)
                          call la_sger(mb,js - 1,one,rhs(3),1,e(1,js),1,f(is,1), &
                                    ldf)
                       end if
                       if (i < p) then
                          call la_sgemv('T',mb,m - ie,-one,a(is,ie + 1),lda,rhs(1),1, &
                                    one,c(ie + 1,js),1)
                          call la_sgemv('T',mb,m - ie,-one,d(is,ie + 1),ldd,rhs(3),1, &
                                    one,c(ie + 1,js),1)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z**t * x = rhs
                       call la_slaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(5,1) = -b(js,js)
                       z(7,1) = -b(jsp1,js)
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(6,2) = -b(js,js)
                       z(8,2) = -b(jsp1,js)
                       z(3,3) = a(is,is)
                       z(4,3) = a(is,isp1)
                       z(5,3) = -b(js,jsp1)
                       z(7,3) = -b(jsp1,jsp1)
                       z(3,4) = a(isp1,is)
                       z(4,4) = a(isp1,isp1)
                       z(6,4) = -b(js,jsp1)
                       z(8,4) = -b(jsp1,jsp1)
                       z(1,5) = d(is,is)
                       z(2,5) = d(is,isp1)
                       z(5,5) = -e(js,js)
                       z(2,6) = d(isp1,isp1)
                       z(6,6) = -e(js,js)
                       z(3,7) = d(is,is)
                       z(4,7) = d(is,isp1)
                       z(5,7) = -e(js,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(4,8) = d(isp1,isp1)
                       z(6,8) = -e(js,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_scopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_scopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z**t * x = rhs
                       call la_sgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_sgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_sscal(m,scaloc,c(1,k),1)
                             call la_sscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_scopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_scopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_sgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1, &
                                    js),ldb,one,f(is,1),ldf)
                          call la_sgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1, &
                                    js),lde,one,f(is,1),ldf)
                       end if
                       if (i < p) then
                          call la_sgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c( &
                                    is,js),ldc,one,c(ie + 1,js),ldc)
                          call la_sgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f( &
                                    is,js),ldf,one,c(ie + 1,js),ldc)
                       end if
                    end if
                 end do loop_190
              end do loop_200
           end if
           return
     end subroutine la_stgsy2
     !> DTGSY2: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                (1)
     !> D * R - L * E = scale * F,
     !> using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)
     !> must be in generalized Schur canonical form, i.e. A, B are upper
     !> quasi triangular and D, E are upper triangular. The solution (R, L)
     !> overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor
     !> chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Z*x = scale*b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**T, Im) ],
     !> Ik is the identity matrix of size k and X**T is the transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> In the process of solving (1), we solve a number of such systems
     !> where Dim(In), Dim(In) = 1 or 2.
     !> If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,
     !> which is equivalent to solve for R and L in
     !> A**T * R  + D**T * L   = scale * C           (3)
     !> R  * B**T + L  * E**T  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> sigma_min(Z) using reverse communication with DLACON.
     !> DTGSY2 also (IJOB >= 1) contributes to the computation in DTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of the matrix pair in
     !> DTGSYL. See DTGSYL for details.

     pure subroutine la_dtgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,iwork,pq,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info,pq
           real(dp),intent(inout) :: rdscal,rdsum
           real(dp),intent(out) :: scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(dp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
        ! replaced various illegal calls to la_dcopy by calls to la_dlaset.
        ! sven hammarling, 27/5/02.
           ! Parameters
           integer(ilp),parameter :: ldz = 8

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ie,ierr,ii,is,isp1,j,je,jj,js,jsp1,k,mb,nb,p,q, &
                     zdim
           real(dp) :: alpha,scaloc
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           real(dp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTGSY2',-info)
              return
           end if
           ! determine block structure of a
           pq = 0
           p = 0
           i = 1
           10 continue
           if (i > m) go to 20
           p = p + 1
           iwork(p) = i
           if (i == m) go to 20
           if (a(i + 1,i) /= zero) then
              i = i + 2
           else
              i = i + 1
           end if
           go to 10
           20 continue
           iwork(p + 1) = m + 1
           ! determine block structure of b
           q = p + 1
           j = 1
           30 continue
           if (j > n) go to 40
           q = q + 1
           iwork(q) = j
           if (j == n) go to 40
           if (b(j + 1,j) /= zero) then
              j = j + 2
           else
              j = j + 1
           end if
           go to 30
           40 continue
           iwork(q + 1) = n + 1
           pq = p*(q - p - 1)
           if (notran) then
              ! solve (i, j) - subsystem
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
              scale = one
              scaloc = one
              loop_120: do j = p + 2,q
                 js = iwork(j)
                 jsp1 = js + 1
                 je = iwork(j + 1) - 1
                 nb = je - js + 1
                 loop_110: do i = p,1,-1
                    is = iwork(i)
                    isp1 = is + 1
                    ie = iwork(i + 1) - 1
                    mb = ie - is + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = d(is,is)
                       z(1,2) = -b(js,js)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_dscal(m,scaloc,c(1,k),1)
                                call la_dscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_dlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          alpha = -rhs(1)
                          call la_daxpy(is - 1,alpha,a(1,is),1,c(1,js),1)
                          call la_daxpy(is - 1,alpha,d(1,is),1,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_daxpy(n - je,rhs(2),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_daxpy(n - je,rhs(2),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = zero
                       z(4,2) = d(is,is)
                       z(1,3) = -b(js,js)
                       z(2,3) = -b(js,jsp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = -e(js,jsp1)
                       z(1,4) = -b(jsp1,js)
                       z(2,4) = -b(jsp1,jsp1)
                       z(3,4) = zero
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_dscal(m,scaloc,c(1,k),1)
                                call la_dscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_dlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_dger(is - 1,nb,-one,a(1,is),1,rhs(1),1,c(1,js), &
                                     ldc)
                          call la_dger(is - 1,nb,-one,d(1,is),1,rhs(1),1,f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_daxpy(n - je,rhs(3),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_daxpy(n - je,rhs(3),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                          call la_daxpy(n - je,rhs(4),b(jsp1,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_daxpy(n - je,rhs(4),e(jsp1,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = d(is,isp1)
                       z(4,2) = d(isp1,isp1)
                       z(1,3) = -b(js,js)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = -b(js,js)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_dscal(m,scaloc,c(1,k),1)
                                call la_dscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_dlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_dgemv('N',is - 1,mb,-one,a(1,is),lda,rhs(1),1, &
                                    one,c(1,js),1)
                          call la_dgemv('N',is - 1,mb,-one,d(1,is),ldd,rhs(1),1, &
                                    one,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_dger(mb,n - je,one,rhs(3),1,b(js,je + 1),ldb,c(is, &
                                    je + 1),ldc)
                          call la_dger(mb,n - je,one,rhs(3),1,e(js,je + 1),lde,f(is, &
                                    je + 1),ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z * x = rhs
                       call la_dlaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(5,1) = d(is,is)
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(5,2) = d(is,isp1)
                       z(6,2) = d(isp1,isp1)
                       z(3,3) = a(is,is)
                       z(4,3) = a(isp1,is)
                       z(7,3) = d(is,is)
                       z(3,4) = a(is,isp1)
                       z(4,4) = a(isp1,isp1)
                       z(7,4) = d(is,isp1)
                       z(8,4) = d(isp1,isp1)
                       z(1,5) = -b(js,js)
                       z(3,5) = -b(js,jsp1)
                       z(5,5) = -e(js,js)
                       z(7,5) = -e(js,jsp1)
                       z(2,6) = -b(js,js)
                       z(4,6) = -b(js,jsp1)
                       z(6,6) = -e(js,js)
                       z(8,6) = -e(js,jsp1)
                       z(1,7) = -b(jsp1,js)
                       z(3,7) = -b(jsp1,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(2,8) = -b(jsp1,js)
                       z(4,8) = -b(jsp1,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_dcopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_dcopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_dscal(m,scaloc,c(1,k),1)
                                call la_dscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_dlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_dcopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_dcopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_dgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,rhs(1 &
                                    ),mb,one,c(1,js),ldc)
                          call la_dgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,rhs(1 &
                                    ),mb,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          k = mb*nb + 1
                          call la_dgemm('N','N',mb,n - je,nb,one,rhs(k),mb,b(js,je + &
                                    1),ldb,one,c(is,je + 1),ldc)
                          call la_dgemm('N','N',mb,n - je,nb,one,rhs(k),mb,e(js,je + &
                                    1),lde,one,f(is,je + 1),ldf)
                       end if
                    end if
                 end do loop_110
              end do loop_120
           else
              ! solve (i, j) - subsystem
                   ! a(i, i)**t * r(i, j) + d(i, i)**t * l(j, j)  =  c(i, j)
                   ! r(i, i)  * b(j, j) + l(i, j)  * e(j, j)  = -f(i, j)
              ! for i = 1, 2, ..., p, j = q, q - 1, ..., 1
              scale = one
              scaloc = one
              loop_200: do i = 1,p
                 is = iwork(i)
                 isp1 = is + 1
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_190: do j = q,p + 2,-1
                    js = iwork(j)
                    jsp1 = js + 1
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = -b(js,js)
                       z(1,2) = d(is,is)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z**t * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_dscal(m,scaloc,c(1,k),1)
                             call la_dscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          alpha = rhs(1)
                          call la_daxpy(js - 1,alpha,b(1,js),1,f(is,1),ldf)
                          alpha = rhs(2)
                          call la_daxpy(js - 1,alpha,e(1,js),1,f(is,1),ldf)
                       end if
                       if (i < p) then
                          alpha = -rhs(1)
                          call la_daxpy(m - ie,alpha,a(is,ie + 1),lda,c(ie + 1,js),1)

                          alpha = -rhs(2)
                          call la_daxpy(m - ie,alpha,d(is,ie + 1),ldd,c(ie + 1,js),1)

                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = -b(js,js)
                       z(4,1) = -b(jsp1,js)
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = -b(js,jsp1)
                       z(4,2) = -b(jsp1,jsp1)
                       z(1,3) = d(is,is)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(is,is)
                       z(3,4) = -e(js,jsp1)
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z**t * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_dscal(m,scaloc,c(1,k),1)
                             call la_dscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_daxpy(js - 1,rhs(1),b(1,js),1,f(is,1),ldf)

                          call la_daxpy(js - 1,rhs(2),b(1,jsp1),1,f(is,1),ldf)

                          call la_daxpy(js - 1,rhs(3),e(1,js),1,f(is,1),ldf)

                          call la_daxpy(js - 1,rhs(4),e(1,jsp1),1,f(is,1),ldf)

                       end if
                       if (i < p) then
                          call la_dger(m - ie,nb,-one,a(is,ie + 1),lda,rhs(1),1,c(ie + &
                                    1,js),ldc)
                          call la_dger(m - ie,nb,-one,d(is,ie + 1),ldd,rhs(3),1,c(ie + &
                                    1,js),ldc)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(3,1) = -b(js,js)
                       z(4,1) = zero
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = zero
                       z(4,2) = -b(js,js)
                       z(1,3) = d(is,is)
                       z(2,3) = d(is,isp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(isp1,isp1)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z**t * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_dscal(m,scaloc,c(1,k),1)
                             call la_dscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_dger(mb,js - 1,one,rhs(1),1,b(1,js),1,f(is,1), &
                                    ldf)
                          call la_dger(mb,js - 1,one,rhs(3),1,e(1,js),1,f(is,1), &
                                    ldf)
                       end if
                       if (i < p) then
                          call la_dgemv('T',mb,m - ie,-one,a(is,ie + 1),lda,rhs(1),1, &
                                    one,c(ie + 1,js),1)
                          call la_dgemv('T',mb,m - ie,-one,d(is,ie + 1),ldd,rhs(3),1, &
                                    one,c(ie + 1,js),1)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z**t * x = rhs
                       call la_dlaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(5,1) = -b(js,js)
                       z(7,1) = -b(jsp1,js)
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(6,2) = -b(js,js)
                       z(8,2) = -b(jsp1,js)
                       z(3,3) = a(is,is)
                       z(4,3) = a(is,isp1)
                       z(5,3) = -b(js,jsp1)
                       z(7,3) = -b(jsp1,jsp1)
                       z(3,4) = a(isp1,is)
                       z(4,4) = a(isp1,isp1)
                       z(6,4) = -b(js,jsp1)
                       z(8,4) = -b(jsp1,jsp1)
                       z(1,5) = d(is,is)
                       z(2,5) = d(is,isp1)
                       z(5,5) = -e(js,js)
                       z(2,6) = d(isp1,isp1)
                       z(6,6) = -e(js,js)
                       z(3,7) = d(is,is)
                       z(4,7) = d(is,isp1)
                       z(5,7) = -e(js,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(4,8) = d(isp1,isp1)
                       z(6,8) = -e(js,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_dcopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_dcopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z**t * x = rhs
                       call la_dgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_dgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_dscal(m,scaloc,c(1,k),1)
                             call la_dscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_dcopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_dcopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_dgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1, &
                                    js),ldb,one,f(is,1),ldf)
                          call la_dgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1, &
                                    js),lde,one,f(is,1),ldf)
                       end if
                       if (i < p) then
                          call la_dgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c( &
                                    is,js),ldc,one,c(ie + 1,js),ldc)
                          call la_dgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f( &
                                    is,js),ldf,one,c(ie + 1,js),ldc)
                       end if
                    end if
                 end do loop_190
              end do loop_200
           end if
           return
     end subroutine la_dtgsy2
#ifdef LA_WITH_XDP
     !> XTGSY2: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                (1)
     !> D * R - L * E = scale * F,
     !> using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)
     !> must be in generalized Schur canonical form, i.e. A, B are upper
     !> quasi triangular and D, E are upper triangular. The solution (R, L)
     !> overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor
     !> chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Z*x = scale*b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**T, Im) ],
     !> Ik is the identity matrix of size k and X**T is the transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> In the process of solving (1), we solve a number of such systems
     !> where Dim(In), Dim(In) = 1 or 2.
     !> If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,
     !> which is equivalent to solve for R and L in
     !> A**T * R  + D**T * L   = scale * C           (3)
     !> R  * B**T + L  * E**T  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> sigma_min(Z) using reverse communication with XLACON.
     !> XTGSY2 also (IJOB >= 1) contributes to the computation in XTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of the matrix pair in
     !> XTGSYL. See XTGSYL for details.

     pure subroutine la_xtgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,iwork,pq,info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info,pq
           real(xdp),intent(inout) :: rdscal,rdsum
           real(xdp),intent(out) :: scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(xdp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
        ! replaced various illegal calls to la_xcopy by calls to la_xlaset.
        ! sven hammarling, 27/5/02.
           ! Parameters
           integer(ilp),parameter :: ldz = 8

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ie,ierr,ii,is,isp1,j,je,jj,js,jsp1,k,mb,nb,p,q, &
                     zdim
           real(xdp) :: alpha,scaloc
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           real(xdp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTGSY2',-info)
              return
           end if
           ! determine block structure of a
           pq = 0
           p = 0
           i = 1
           10 continue
           if (i > m) go to 20
           p = p + 1
           iwork(p) = i
           if (i == m) go to 20
           if (a(i + 1,i) /= zero) then
              i = i + 2
           else
              i = i + 1
           end if
           go to 10
           20 continue
           iwork(p + 1) = m + 1
           ! determine block structure of b
           q = p + 1
           j = 1
           30 continue
           if (j > n) go to 40
           q = q + 1
           iwork(q) = j
           if (j == n) go to 40
           if (b(j + 1,j) /= zero) then
              j = j + 2
           else
              j = j + 1
           end if
           go to 30
           40 continue
           iwork(q + 1) = n + 1
           pq = p*(q - p - 1)
           if (notran) then
              ! solve (i, j) - subsystem
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
              scale = one
              scaloc = one
              loop_120: do j = p + 2,q
                 js = iwork(j)
                 jsp1 = js + 1
                 je = iwork(j + 1) - 1
                 nb = je - js + 1
                 loop_110: do i = p,1,-1
                    is = iwork(i)
                    isp1 = is + 1
                    ie = iwork(i + 1) - 1
                    mb = ie - is + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = d(is,is)
                       z(1,2) = -b(js,js)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_xscal(m,scaloc,c(1,k),1)
                                call la_xscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_xlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          alpha = -rhs(1)
                          call la_xaxpy(is - 1,alpha,a(1,is),1,c(1,js),1)
                          call la_xaxpy(is - 1,alpha,d(1,is),1,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_xaxpy(n - je,rhs(2),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_xaxpy(n - je,rhs(2),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = zero
                       z(4,2) = d(is,is)
                       z(1,3) = -b(js,js)
                       z(2,3) = -b(js,jsp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = -e(js,jsp1)
                       z(1,4) = -b(jsp1,js)
                       z(2,4) = -b(jsp1,jsp1)
                       z(3,4) = zero
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_xscal(m,scaloc,c(1,k),1)
                                call la_xscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_xlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_xger(is - 1,nb,-one,a(1,is),1,rhs(1),1,c(1,js), &
                                     ldc)
                          call la_xger(is - 1,nb,-one,d(1,is),1,rhs(1),1,f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_xaxpy(n - je,rhs(3),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_xaxpy(n - je,rhs(3),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                          call la_xaxpy(n - je,rhs(4),b(jsp1,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_xaxpy(n - je,rhs(4),e(jsp1,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = d(is,isp1)
                       z(4,2) = d(isp1,isp1)
                       z(1,3) = -b(js,js)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = -b(js,js)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_xscal(m,scaloc,c(1,k),1)
                                call la_xscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_xlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_xgemv('N',is - 1,mb,-one,a(1,is),lda,rhs(1),1, &
                                    one,c(1,js),1)
                          call la_xgemv('N',is - 1,mb,-one,d(1,is),ldd,rhs(1),1, &
                                    one,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_xger(mb,n - je,one,rhs(3),1,b(js,je + 1),ldb,c(is, &
                                    je + 1),ldc)
                          call la_xger(mb,n - je,one,rhs(3),1,e(js,je + 1),lde,f(is, &
                                    je + 1),ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z * x = rhs
                       call la_xlaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(5,1) = d(is,is)
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(5,2) = d(is,isp1)
                       z(6,2) = d(isp1,isp1)
                       z(3,3) = a(is,is)
                       z(4,3) = a(isp1,is)
                       z(7,3) = d(is,is)
                       z(3,4) = a(is,isp1)
                       z(4,4) = a(isp1,isp1)
                       z(7,4) = d(is,isp1)
                       z(8,4) = d(isp1,isp1)
                       z(1,5) = -b(js,js)
                       z(3,5) = -b(js,jsp1)
                       z(5,5) = -e(js,js)
                       z(7,5) = -e(js,jsp1)
                       z(2,6) = -b(js,js)
                       z(4,6) = -b(js,jsp1)
                       z(6,6) = -e(js,js)
                       z(8,6) = -e(js,jsp1)
                       z(1,7) = -b(jsp1,js)
                       z(3,7) = -b(jsp1,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(2,8) = -b(jsp1,js)
                       z(4,8) = -b(jsp1,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_xcopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_xcopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_xscal(m,scaloc,c(1,k),1)
                                call la_xscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_xlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_xcopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_xcopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_xgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,rhs(1 &
                                    ),mb,one,c(1,js),ldc)
                          call la_xgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,rhs(1 &
                                    ),mb,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          k = mb*nb + 1
                          call la_xgemm('N','N',mb,n - je,nb,one,rhs(k),mb,b(js,je + &
                                    1),ldb,one,c(is,je + 1),ldc)
                          call la_xgemm('N','N',mb,n - je,nb,one,rhs(k),mb,e(js,je + &
                                    1),lde,one,f(is,je + 1),ldf)
                       end if
                    end if
                 end do loop_110
              end do loop_120
           else
              ! solve (i, j) - subsystem
                   ! a(i, i)**t * r(i, j) + d(i, i)**t * l(j, j)  =  c(i, j)
                   ! r(i, i)  * b(j, j) + l(i, j)  * e(j, j)  = -f(i, j)
              ! for i = 1, 2, ..., p, j = q, q - 1, ..., 1
              scale = one
              scaloc = one
              loop_200: do i = 1,p
                 is = iwork(i)
                 isp1 = is + 1
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_190: do j = q,p + 2,-1
                    js = iwork(j)
                    jsp1 = js + 1
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = -b(js,js)
                       z(1,2) = d(is,is)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z**t * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_xscal(m,scaloc,c(1,k),1)
                             call la_xscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          alpha = rhs(1)
                          call la_xaxpy(js - 1,alpha,b(1,js),1,f(is,1),ldf)
                          alpha = rhs(2)
                          call la_xaxpy(js - 1,alpha,e(1,js),1,f(is,1),ldf)
                       end if
                       if (i < p) then
                          alpha = -rhs(1)
                          call la_xaxpy(m - ie,alpha,a(is,ie + 1),lda,c(ie + 1,js),1)

                          alpha = -rhs(2)
                          call la_xaxpy(m - ie,alpha,d(is,ie + 1),ldd,c(ie + 1,js),1)

                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = -b(js,js)
                       z(4,1) = -b(jsp1,js)
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = -b(js,jsp1)
                       z(4,2) = -b(jsp1,jsp1)
                       z(1,3) = d(is,is)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(is,is)
                       z(3,4) = -e(js,jsp1)
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z**t * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_xscal(m,scaloc,c(1,k),1)
                             call la_xscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_xaxpy(js - 1,rhs(1),b(1,js),1,f(is,1),ldf)

                          call la_xaxpy(js - 1,rhs(2),b(1,jsp1),1,f(is,1),ldf)

                          call la_xaxpy(js - 1,rhs(3),e(1,js),1,f(is,1),ldf)

                          call la_xaxpy(js - 1,rhs(4),e(1,jsp1),1,f(is,1),ldf)

                       end if
                       if (i < p) then
                          call la_xger(m - ie,nb,-one,a(is,ie + 1),lda,rhs(1),1,c(ie + &
                                    1,js),ldc)
                          call la_xger(m - ie,nb,-one,d(is,ie + 1),ldd,rhs(3),1,c(ie + &
                                    1,js),ldc)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(3,1) = -b(js,js)
                       z(4,1) = zero
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = zero
                       z(4,2) = -b(js,js)
                       z(1,3) = d(is,is)
                       z(2,3) = d(is,isp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(isp1,isp1)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z**t * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_xscal(m,scaloc,c(1,k),1)
                             call la_xscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_xger(mb,js - 1,one,rhs(1),1,b(1,js),1,f(is,1), &
                                    ldf)
                          call la_xger(mb,js - 1,one,rhs(3),1,e(1,js),1,f(is,1), &
                                    ldf)
                       end if
                       if (i < p) then
                          call la_xgemv('T',mb,m - ie,-one,a(is,ie + 1),lda,rhs(1),1, &
                                    one,c(ie + 1,js),1)
                          call la_xgemv('T',mb,m - ie,-one,d(is,ie + 1),ldd,rhs(3),1, &
                                    one,c(ie + 1,js),1)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z**t * x = rhs
                       call la_xlaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(5,1) = -b(js,js)
                       z(7,1) = -b(jsp1,js)
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(6,2) = -b(js,js)
                       z(8,2) = -b(jsp1,js)
                       z(3,3) = a(is,is)
                       z(4,3) = a(is,isp1)
                       z(5,3) = -b(js,jsp1)
                       z(7,3) = -b(jsp1,jsp1)
                       z(3,4) = a(isp1,is)
                       z(4,4) = a(isp1,isp1)
                       z(6,4) = -b(js,jsp1)
                       z(8,4) = -b(jsp1,jsp1)
                       z(1,5) = d(is,is)
                       z(2,5) = d(is,isp1)
                       z(5,5) = -e(js,js)
                       z(2,6) = d(isp1,isp1)
                       z(6,6) = -e(js,js)
                       z(3,7) = d(is,is)
                       z(4,7) = d(is,isp1)
                       z(5,7) = -e(js,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(4,8) = d(isp1,isp1)
                       z(6,8) = -e(js,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_xcopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_xcopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z**t * x = rhs
                       call la_xgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_xgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_xscal(m,scaloc,c(1,k),1)
                             call la_xscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_xcopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_xcopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_xgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1, &
                                    js),ldb,one,f(is,1),ldf)
                          call la_xgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1, &
                                    js),lde,one,f(is,1),ldf)
                       end if
                       if (i < p) then
                          call la_xgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c( &
                                    is,js),ldc,one,c(ie + 1,js),ldc)
                          call la_xgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f( &
                                    is,js),ldf,one,c(ie + 1,js),ldc)
                       end if
                    end if
                 end do loop_190
              end do loop_200
           end if
           return
     end subroutine la_xtgsy2
#endif
#ifdef LA_WITH_QP
     !> QTGSY2: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                (1)
     !> D * R - L * E = scale * F,
     !> using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)
     !> must be in generalized Schur canonical form, i.e. A, B are upper
     !> quasi triangular and D, E are upper triangular. The solution (R, L)
     !> overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor
     !> chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Z*x = scale*b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**T, Im) ],
     !> Ik is the identity matrix of size k and X**T is the transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> In the process of solving (1), we solve a number of such systems
     !> where Dim(In), Dim(In) = 1 or 2.
     !> If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,
     !> which is equivalent to solve for R and L in
     !> A**T * R  + D**T * L   = scale * C           (3)
     !> R  * B**T + L  * E**T  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> sigma_min(Z) using reverse communication with QLACON.
     !> QTGSY2 also (IJOB >= 1) contributes to the computation in QTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of the matrix pair in
     !> QTGSYL. See QTGSYL for details.

     pure subroutine la_qtgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,iwork,pq,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info,pq
           real(qp),intent(inout) :: rdscal,rdsum
           real(qp),intent(out) :: scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(qp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
        ! replaced various illegal calls to la_qcopy by calls to la_qlaset.
        ! sven hammarling, 27/5/02.
           ! Parameters
           integer(ilp),parameter :: ldz = 8

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ie,ierr,ii,is,isp1,j,je,jj,js,jsp1,k,mb,nb,p,q, &
                     zdim
           real(qp) :: alpha,scaloc
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           real(qp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTGSY2',-info)
              return
           end if
           ! determine block structure of a
           pq = 0
           p = 0
           i = 1
           10 continue
           if (i > m) go to 20
           p = p + 1
           iwork(p) = i
           if (i == m) go to 20
           if (a(i + 1,i) /= zero) then
              i = i + 2
           else
              i = i + 1
           end if
           go to 10
           20 continue
           iwork(p + 1) = m + 1
           ! determine block structure of b
           q = p + 1
           j = 1
           30 continue
           if (j > n) go to 40
           q = q + 1
           iwork(q) = j
           if (j == n) go to 40
           if (b(j + 1,j) /= zero) then
              j = j + 2
           else
              j = j + 1
           end if
           go to 30
           40 continue
           iwork(q + 1) = n + 1
           pq = p*(q - p - 1)
           if (notran) then
              ! solve (i, j) - subsystem
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
              scale = one
              scaloc = one
              loop_120: do j = p + 2,q
                 js = iwork(j)
                 jsp1 = js + 1
                 je = iwork(j + 1) - 1
                 nb = je - js + 1
                 loop_110: do i = p,1,-1
                    is = iwork(i)
                    isp1 = is + 1
                    ie = iwork(i + 1) - 1
                    mb = ie - is + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = d(is,is)
                       z(1,2) = -b(js,js)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_qscal(m,scaloc,c(1,k),1)
                                call la_qscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_qlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          alpha = -rhs(1)
                          call la_qaxpy(is - 1,alpha,a(1,is),1,c(1,js),1)
                          call la_qaxpy(is - 1,alpha,d(1,is),1,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_qaxpy(n - je,rhs(2),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_qaxpy(n - je,rhs(2),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = zero
                       z(4,2) = d(is,is)
                       z(1,3) = -b(js,js)
                       z(2,3) = -b(js,jsp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = -e(js,jsp1)
                       z(1,4) = -b(jsp1,js)
                       z(2,4) = -b(jsp1,jsp1)
                       z(3,4) = zero
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_qscal(m,scaloc,c(1,k),1)
                                call la_qscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_qlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_qger(is - 1,nb,-one,a(1,is),1,rhs(1),1,c(1,js), &
                                     ldc)
                          call la_qger(is - 1,nb,-one,d(1,is),1,rhs(1),1,f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_qaxpy(n - je,rhs(3),b(js,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_qaxpy(n - je,rhs(3),e(js,je + 1),lde,f(is,je + 1), &
                                    ldf)
                          call la_qaxpy(n - je,rhs(4),b(jsp1,je + 1),ldb,c(is,je + 1), &
                                    ldc)
                          call la_qaxpy(n - je,rhs(4),e(jsp1,je + 1),lde,f(is,je + 1), &
                                    ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(3,1) = d(is,is)
                       z(4,1) = zero
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = d(is,isp1)
                       z(4,2) = d(isp1,isp1)
                       z(1,3) = -b(js,js)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = -b(js,js)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_qscal(m,scaloc,c(1,k),1)
                                call la_qscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_qlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_qgemv('N',is - 1,mb,-one,a(1,is),lda,rhs(1),1, &
                                    one,c(1,js),1)
                          call la_qgemv('N',is - 1,mb,-one,d(1,is),ldd,rhs(1),1, &
                                    one,f(1,js),1)
                       end if
                       if (j < q) then
                          call la_qger(mb,n - je,one,rhs(3),1,b(js,je + 1),ldb,c(is, &
                                    je + 1),ldc)
                          call la_qger(mb,n - je,one,rhs(3),1,e(js,je + 1),lde,f(is, &
                                    je + 1),ldf)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z * x = rhs
                       call la_qlaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(isp1,is)
                       z(5,1) = d(is,is)
                       z(1,2) = a(is,isp1)
                       z(2,2) = a(isp1,isp1)
                       z(5,2) = d(is,isp1)
                       z(6,2) = d(isp1,isp1)
                       z(3,3) = a(is,is)
                       z(4,3) = a(isp1,is)
                       z(7,3) = d(is,is)
                       z(3,4) = a(is,isp1)
                       z(4,4) = a(isp1,isp1)
                       z(7,4) = d(is,isp1)
                       z(8,4) = d(isp1,isp1)
                       z(1,5) = -b(js,js)
                       z(3,5) = -b(js,jsp1)
                       z(5,5) = -e(js,js)
                       z(7,5) = -e(js,jsp1)
                       z(2,6) = -b(js,js)
                       z(4,6) = -b(js,jsp1)
                       z(6,6) = -e(js,js)
                       z(8,6) = -e(js,jsp1)
                       z(1,7) = -b(jsp1,js)
                       z(3,7) = -b(jsp1,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(2,8) = -b(jsp1,js)
                       z(4,8) = -b(jsp1,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_qcopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_qcopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       if (ijob == 0) then
                          call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                          if (scaloc /= one) then
                             do k = 1,n
                                call la_qscal(m,scaloc,c(1,k),1)
                                call la_qscal(m,scaloc,f(1,k),1)
                             end do
                             scale = scale*scaloc
                          end if
                       else
                          call la_qlatdf(ijob,zdim,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_qcopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_qcopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_qgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,rhs(1 &
                                    ),mb,one,c(1,js),ldc)
                          call la_qgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,rhs(1 &
                                    ),mb,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          k = mb*nb + 1
                          call la_qgemm('N','N',mb,n - je,nb,one,rhs(k),mb,b(js,je + &
                                    1),ldb,one,c(is,je + 1),ldc)
                          call la_qgemm('N','N',mb,n - je,nb,one,rhs(k),mb,e(js,je + &
                                    1),lde,one,f(is,je + 1),ldf)
                       end if
                    end if
                 end do loop_110
              end do loop_120
           else
              ! solve (i, j) - subsystem
                   ! a(i, i)**t * r(i, j) + d(i, i)**t * l(j, j)  =  c(i, j)
                   ! r(i, i)  * b(j, j) + l(i, j)  * e(j, j)  = -f(i, j)
              ! for i = 1, 2, ..., p, j = q, q - 1, ..., 1
              scale = one
              scaloc = one
              loop_200: do i = 1,p
                 is = iwork(i)
                 isp1 = is + 1
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_190: do j = q,p + 2,-1
                    js = iwork(j)
                    jsp1 = js + 1
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    zdim = mb*nb*2
                    if ((mb == 1) .and. (nb == 1)) then
                       ! build a 2-by-2 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = -b(js,js)
                       z(1,2) = d(is,is)
                       z(2,2) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = f(is,js)
                       ! solve z**t * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_qscal(m,scaloc,c(1,k),1)
                             call la_qscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       f(is,js) = rhs(2)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          alpha = rhs(1)
                          call la_qaxpy(js - 1,alpha,b(1,js),1,f(is,1),ldf)
                          alpha = rhs(2)
                          call la_qaxpy(js - 1,alpha,e(1,js),1,f(is,1),ldf)
                       end if
                       if (i < p) then
                          alpha = -rhs(1)
                          call la_qaxpy(m - ie,alpha,a(is,ie + 1),lda,c(ie + 1,js),1)

                          alpha = -rhs(2)
                          call la_qaxpy(m - ie,alpha,d(is,ie + 1),ldd,c(ie + 1,js),1)

                       end if
                    else if ((mb == 1) .and. (nb == 2)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = zero
                       z(3,1) = -b(js,js)
                       z(4,1) = -b(jsp1,js)
                       z(1,2) = zero
                       z(2,2) = a(is,is)
                       z(3,2) = -b(js,jsp1)
                       z(4,2) = -b(jsp1,jsp1)
                       z(1,3) = d(is,is)
                       z(2,3) = zero
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(is,is)
                       z(3,4) = -e(js,jsp1)
                       z(4,4) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(is,jsp1)
                       rhs(3) = f(is,js)
                       rhs(4) = f(is,jsp1)
                       ! solve z**t * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_qscal(m,scaloc,c(1,k),1)
                             call la_qscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(is,jsp1) = rhs(2)
                       f(is,js) = rhs(3)
                       f(is,jsp1) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_qaxpy(js - 1,rhs(1),b(1,js),1,f(is,1),ldf)

                          call la_qaxpy(js - 1,rhs(2),b(1,jsp1),1,f(is,1),ldf)

                          call la_qaxpy(js - 1,rhs(3),e(1,js),1,f(is,1),ldf)

                          call la_qaxpy(js - 1,rhs(4),e(1,jsp1),1,f(is,1),ldf)

                       end if
                       if (i < p) then
                          call la_qger(m - ie,nb,-one,a(is,ie + 1),lda,rhs(1),1,c(ie + &
                                    1,js),ldc)
                          call la_qger(m - ie,nb,-one,d(is,ie + 1),ldd,rhs(3),1,c(ie + &
                                    1,js),ldc)
                       end if
                    else if ((mb == 2) .and. (nb == 1)) then
                       ! build a 4-by-4 system z**t * x = rhs
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(3,1) = -b(js,js)
                       z(4,1) = zero
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(3,2) = zero
                       z(4,2) = -b(js,js)
                       z(1,3) = d(is,is)
                       z(2,3) = d(is,isp1)
                       z(3,3) = -e(js,js)
                       z(4,3) = zero
                       z(1,4) = zero
                       z(2,4) = d(isp1,isp1)
                       z(3,4) = zero
                       z(4,4) = -e(js,js)
                       ! set up right hand side(s)
                       rhs(1) = c(is,js)
                       rhs(2) = c(isp1,js)
                       rhs(3) = f(is,js)
                       rhs(4) = f(isp1,js)
                       ! solve z**t * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_qscal(m,scaloc,c(1,k),1)
                             call la_qscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       c(is,js) = rhs(1)
                       c(isp1,js) = rhs(2)
                       f(is,js) = rhs(3)
                       f(isp1,js) = rhs(4)
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_qger(mb,js - 1,one,rhs(1),1,b(1,js),1,f(is,1), &
                                    ldf)
                          call la_qger(mb,js - 1,one,rhs(3),1,e(1,js),1,f(is,1), &
                                    ldf)
                       end if
                       if (i < p) then
                          call la_qgemv('T',mb,m - ie,-one,a(is,ie + 1),lda,rhs(1),1, &
                                    one,c(ie + 1,js),1)
                          call la_qgemv('T',mb,m - ie,-one,d(is,ie + 1),ldd,rhs(3),1, &
                                    one,c(ie + 1,js),1)
                       end if
                    else if ((mb == 2) .and. (nb == 2)) then
                       ! build an 8-by-8 system z**t * x = rhs
                       call la_qlaset('F',ldz,ldz,zero,zero,z,ldz)
                       z(1,1) = a(is,is)
                       z(2,1) = a(is,isp1)
                       z(5,1) = -b(js,js)
                       z(7,1) = -b(jsp1,js)
                       z(1,2) = a(isp1,is)
                       z(2,2) = a(isp1,isp1)
                       z(6,2) = -b(js,js)
                       z(8,2) = -b(jsp1,js)
                       z(3,3) = a(is,is)
                       z(4,3) = a(is,isp1)
                       z(5,3) = -b(js,jsp1)
                       z(7,3) = -b(jsp1,jsp1)
                       z(3,4) = a(isp1,is)
                       z(4,4) = a(isp1,isp1)
                       z(6,4) = -b(js,jsp1)
                       z(8,4) = -b(jsp1,jsp1)
                       z(1,5) = d(is,is)
                       z(2,5) = d(is,isp1)
                       z(5,5) = -e(js,js)
                       z(2,6) = d(isp1,isp1)
                       z(6,6) = -e(js,js)
                       z(3,7) = d(is,is)
                       z(4,7) = d(is,isp1)
                       z(5,7) = -e(js,jsp1)
                       z(7,7) = -e(jsp1,jsp1)
                       z(4,8) = d(isp1,isp1)
                       z(6,8) = -e(js,jsp1)
                       z(8,8) = -e(jsp1,jsp1)
                       ! set up right hand side(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_qcopy(mb,c(is,js + jj),1,rhs(k),1)
                          call la_qcopy(mb,f(is,js + jj),1,rhs(ii),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! solve z**t * x = rhs
                       call la_qgetc2(zdim,z,ldz,ipiv,jpiv,ierr)
                       if (ierr > 0) info = ierr
                       call la_qgesc2(zdim,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_qscal(m,scaloc,c(1,k),1)
                             call la_qscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! unpack solution vector(s)
                       k = 1
                       ii = mb*nb + 1
                       do jj = 0,nb - 1
                          call la_qcopy(mb,rhs(k),1,c(is,js + jj),1)
                          call la_qcopy(mb,rhs(ii),1,f(is,js + jj),1)
                          k = k + mb
                          ii = ii + mb
                       end do
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (j > p + 2) then
                          call la_qgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1, &
                                    js),ldb,one,f(is,1),ldf)
                          call la_qgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1, &
                                    js),lde,one,f(is,1),ldf)
                       end if
                       if (i < p) then
                          call la_qgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c( &
                                    is,js),ldc,one,c(ie + 1,js),ldc)
                          call la_qgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f( &
                                    is,js),ldf,one,c(ie + 1,js),ldc)
                       end if
                    end if
                 end do loop_190
              end do loop_200
           end if
           return
     end subroutine la_qtgsy2
#endif

     !> STGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                 (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with real entries. (A, D) and (B, E) must be in
     !> generalized (real) Schur canonical form, i.e. A, B are upper quasi
     !> triangular and D, E are upper triangular.
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve  Zx = scale b, where
     !> Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
     !> [ kron(In, D)  -kron(E**T, Im) ].
     !> Here Ik is the identity matrix of size k and X**T is the transpose of
     !> X. kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'T', STGSYL solves the transposed system Z**T*y = scale*b,
     !> which is equivalent to solve for R and L in
     !> A**T * R + D**T * L = scale * C           (3)
     !> R * B**T + L * E**T = scale * -F
     !> This case (TRANS = 'T') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using SLACON.
     !> If IJOB >= 1, STGSYL computes a Frobenius norm-based estimate
     !> of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z. See [1-2] for more
     !> information.
     !> This is a level 3 BLAS algorithm.

     pure subroutine la_stgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(sp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(sp),intent(inout) :: c(ldc,*),f(ldf,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_scopy by calls to la_slaset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,ppqq,pq,q
           real(sp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: max,real,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine optimal block sizes mb and nb
           mb = la_ilaenv(2,'STGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'STGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_slaset('F',m,n,zero,zero,c,ldc)
                 call la_slaset('F',m,n,zero,zero,f,ldf)
              else if (ijob >= 1 .and. notran) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              loop_30: do iround = 1,isolve
                 ! use unblocked level 2 solver
                 dscale = zero
                 dsum = one
                 pq = 0
                 call la_stgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,iwork,pq,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=sp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=sp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_slacpy('F',m,n,c,ldc,work,m)
                    call la_slacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_slaset('F',m,n,zero,zero,c,ldc)
                    call la_slaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_slacpy('F',m,n,work,m,c,ldc)
                    call la_slacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           if (a(i,i - 1) /= zero) i = i + 1
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           if (b(j,j - 1) /= zero) j = j + 1
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j)-subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1,..., 1; j = 1, 2,..., q
                 dscale = zero
                 dsum = one
                 pq = 0
                 scale = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       ppqq = 0
                       call la_stgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,iwork(q + 2),ppqq,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + ppqq
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_sscal(m,scaloc,c(1,k),1)
                             call la_sscal(m,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_sscal(is - 1,scaloc,c(1,k),1)
                             call la_sscal(is - 1,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_sscal(m - ie,scaloc,c(ie + 1,k),1)
                             call la_sscal(m - ie,scaloc,f(ie + 1,k),1)
                          end do
                          do k = je + 1,n
                             call la_sscal(m,scaloc,c(1,k),1)
                             call la_sscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_sgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,c(is, &
                                    js),ldc,one,c(1,js),ldc)
                          call la_sgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,c(is, &
                                    js),ldc,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          call la_sgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,b(js, &
                                    je + 1),ldb,one,c(is,je + 1),ldc)
                          call la_sgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,e(js, &
                                    je + 1),lde,one,f(is,je + 1),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=sp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=sp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_slacpy('F',m,n,c,ldc,work,m)
                    call la_slacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_slaset('F',m,n,zero,zero,c,ldc)
                    call la_slaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_slacpy('F',m,n,work,m,c,ldc)
                    call la_slacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                   ! a(i, i)**t * r(i, j)  + d(i, i)**t * l(i, j)  =  c(i, j)
                   ! r(i, j)  * b(j, j)**t + l(i, j)  * e(j, j)**t = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_stgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,iwork(q + 2),ppqq,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_sscal(m,scaloc,c(1,k),1)
                          call la_sscal(m,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_sscal(is - 1,scaloc,c(1,k),1)
                          call la_sscal(is - 1,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_sscal(m - ie,scaloc,c(ie + 1,k),1)
                          call la_sscal(m - ie,scaloc,f(ie + 1,k),1)
                       end do
                       do k = je + 1,n
                          call la_sscal(m,scaloc,c(1,k),1)
                          call la_sscal(m,scaloc,f(1,k),1)
                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (j > p + 2) then
                       call la_sgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1,js) &
                                 ,ldb,one,f(is,1),ldf)
                       call la_sgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1,js) &
                                 ,lde,one,f(is,1),ldf)
                    end if
                    if (i < p) then
                       call la_sgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c(is, &
                                 js),ldc,one,c(ie + 1,js),ldc)
                       call la_sgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f(is, &
                                 js),ldf,one,c(ie + 1,js),ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_stgsyl
     !> DTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                 (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with real entries. (A, D) and (B, E) must be in
     !> generalized (real) Schur canonical form, i.e. A, B are upper quasi
     !> triangular and D, E are upper triangular.
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve  Zx = scale b, where
     !> Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
     !> [ kron(In, D)  -kron(E**T, Im) ].
     !> Here Ik is the identity matrix of size k and X**T is the transpose of
     !> X. kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'T', DTGSYL solves the transposed system Z**T*y = scale*b,
     !> which is equivalent to solve for R and L in
     !> A**T * R + D**T * L = scale * C           (3)
     !> R * B**T + L * E**T = scale * -F
     !> This case (TRANS = 'T') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using DLACON.
     !> If IJOB >= 1, DTGSYL computes a Frobenius norm-based estimate
     !> of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z. See [1-2] for more
     !> information.
     !> This is a level 3 BLAS algorithm.

     pure subroutine la_dtgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(dp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(dp),intent(inout) :: c(ldc,*),f(ldf,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_dcopy by calls to la_dlaset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,ppqq,pq,q
           real(dp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: real,max,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine optimal block sizes mb and nb
           mb = la_ilaenv(2,'DTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'DTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_dlaset('F',m,n,zero,zero,c,ldc)
                 call la_dlaset('F',m,n,zero,zero,f,ldf)
              else if (ijob >= 1) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              loop_30: do iround = 1,isolve
                 ! use unblocked level 2 solver
                 dscale = zero
                 dsum = one
                 pq = 0
                 call la_dtgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,iwork,pq,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=dp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=dp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_dlacpy('F',m,n,c,ldc,work,m)
                    call la_dlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_dlaset('F',m,n,zero,zero,c,ldc)
                    call la_dlaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_dlacpy('F',m,n,work,m,c,ldc)
                    call la_dlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           if (a(i,i - 1) /= zero) i = i + 1
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           if (b(j,j - 1) /= zero) j = j + 1
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j)-subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1,..., 1; j = 1, 2,..., q
                 dscale = zero
                 dsum = one
                 pq = 0
                 scale = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       ppqq = 0
                       call la_dtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,iwork(q + 2),ppqq,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + ppqq
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_dscal(m,scaloc,c(1,k),1)
                             call la_dscal(m,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_dscal(is - 1,scaloc,c(1,k),1)
                             call la_dscal(is - 1,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_dscal(m - ie,scaloc,c(ie + 1,k),1)
                             call la_dscal(m - ie,scaloc,f(ie + 1,k),1)
                          end do
                          do k = je + 1,n
                             call la_dscal(m,scaloc,c(1,k),1)
                             call la_dscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_dgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,c(is, &
                                    js),ldc,one,c(1,js),ldc)
                          call la_dgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,c(is, &
                                    js),ldc,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          call la_dgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,b(js, &
                                    je + 1),ldb,one,c(is,je + 1),ldc)
                          call la_dgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,e(js, &
                                    je + 1),lde,one,f(is,je + 1),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=dp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=dp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_dlacpy('F',m,n,c,ldc,work,m)
                    call la_dlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_dlaset('F',m,n,zero,zero,c,ldc)
                    call la_dlaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_dlacpy('F',m,n,work,m,c,ldc)
                    call la_dlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                   ! a(i, i)**t * r(i, j)  + d(i, i)**t * l(i, j)  =  c(i, j)
                   ! r(i, j)  * b(j, j)**t + l(i, j)  * e(j, j)**t = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_dtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,iwork(q + 2),ppqq,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_dscal(m,scaloc,c(1,k),1)
                          call la_dscal(m,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_dscal(is - 1,scaloc,c(1,k),1)
                          call la_dscal(is - 1,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_dscal(m - ie,scaloc,c(ie + 1,k),1)
                          call la_dscal(m - ie,scaloc,f(ie + 1,k),1)
                       end do
                       do k = je + 1,n
                          call la_dscal(m,scaloc,c(1,k),1)
                          call la_dscal(m,scaloc,f(1,k),1)
                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (j > p + 2) then
                       call la_dgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1,js) &
                                 ,ldb,one,f(is,1),ldf)
                       call la_dgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1,js) &
                                 ,lde,one,f(is,1),ldf)
                    end if
                    if (i < p) then
                       call la_dgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c(is, &
                                 js),ldc,one,c(ie + 1,js),ldc)
                       call la_dgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f(is, &
                                 js),ldf,one,c(ie + 1,js),ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_dtgsyl
#ifdef LA_WITH_XDP
     !> XTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                 (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with real entries. (A, D) and (B, E) must be in
     !> generalized (real) Schur canonical form, i.e. A, B are upper quasi
     !> triangular and D, E are upper triangular.
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve  Zx = scale b, where
     !> Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
     !> [ kron(In, D)  -kron(E**T, Im) ].
     !> Here Ik is the identity matrix of size k and X**T is the transpose of
     !> X. kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'T', XTGSYL solves the transposed system Z**T*y = scale*b,
     !> which is equivalent to solve for R and L in
     !> A**T * R + D**T * L = scale * C           (3)
     !> R * B**T + L * E**T = scale * -F
     !> This case (TRANS = 'T') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using XLACON.
     !> If IJOB >= 1, XTGSYL computes a Frobenius norm-based estimate
     !> of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z. See [1-2] for more
     !> information.
     !> This is a level 3 BLAS algorithm.

     pure subroutine la_xtgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(xdp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(xdp),intent(inout) :: c(ldc,*),f(ldf,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_xcopy by calls to la_xlaset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,ppqq,pq,q
           real(xdp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: real,max,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine optimal block sizes mb and nb
           mb = la_ilaenv(2,'XTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'XTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_xlaset('F',m,n,zero,zero,c,ldc)
                 call la_xlaset('F',m,n,zero,zero,f,ldf)
              else if (ijob >= 1) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              loop_30: do iround = 1,isolve
                 ! use unblocked level 2 solver
                 dscale = zero
                 dsum = one
                 pq = 0
                 call la_xtgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,iwork,pq,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=xdp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=xdp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_xlacpy('F',m,n,c,ldc,work,m)
                    call la_xlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_xlaset('F',m,n,zero,zero,c,ldc)
                    call la_xlaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_xlacpy('F',m,n,work,m,c,ldc)
                    call la_xlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           if (a(i,i - 1) /= zero) i = i + 1
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           if (b(j,j - 1) /= zero) j = j + 1
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j)-subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1,..., 1; j = 1, 2,..., q
                 dscale = zero
                 dsum = one
                 pq = 0
                 scale = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       ppqq = 0
                       call la_xtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,iwork(q + 2),ppqq,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + ppqq
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_xscal(m,scaloc,c(1,k),1)
                             call la_xscal(m,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_xscal(is - 1,scaloc,c(1,k),1)
                             call la_xscal(is - 1,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_xscal(m - ie,scaloc,c(ie + 1,k),1)
                             call la_xscal(m - ie,scaloc,f(ie + 1,k),1)
                          end do
                          do k = je + 1,n
                             call la_xscal(m,scaloc,c(1,k),1)
                             call la_xscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_xgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,c(is, &
                                    js),ldc,one,c(1,js),ldc)
                          call la_xgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,c(is, &
                                    js),ldc,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          call la_xgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,b(js, &
                                    je + 1),ldb,one,c(is,je + 1),ldc)
                          call la_xgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,e(js, &
                                    je + 1),lde,one,f(is,je + 1),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=xdp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=xdp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_xlacpy('F',m,n,c,ldc,work,m)
                    call la_xlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_xlaset('F',m,n,zero,zero,c,ldc)
                    call la_xlaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_xlacpy('F',m,n,work,m,c,ldc)
                    call la_xlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                   ! a(i, i)**t * r(i, j)  + d(i, i)**t * l(i, j)  =  c(i, j)
                   ! r(i, j)  * b(j, j)**t + l(i, j)  * e(j, j)**t = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_xtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,iwork(q + 2),ppqq,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_xscal(m,scaloc,c(1,k),1)
                          call la_xscal(m,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_xscal(is - 1,scaloc,c(1,k),1)
                          call la_xscal(is - 1,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_xscal(m - ie,scaloc,c(ie + 1,k),1)
                          call la_xscal(m - ie,scaloc,f(ie + 1,k),1)
                       end do
                       do k = je + 1,n
                          call la_xscal(m,scaloc,c(1,k),1)
                          call la_xscal(m,scaloc,f(1,k),1)
                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (j > p + 2) then
                       call la_xgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1,js) &
                                 ,ldb,one,f(is,1),ldf)
                       call la_xgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1,js) &
                                 ,lde,one,f(is,1),ldf)
                    end if
                    if (i < p) then
                       call la_xgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c(is, &
                                 js),ldc,one,c(ie + 1,js),ldc)
                       call la_xgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f(is, &
                                 js),ldf,one,c(ie + 1,js),ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_xtgsyl
#endif
#ifdef LA_WITH_QP
     !> QTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C                 (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with real entries. (A, D) and (B, E) must be in
     !> generalized (real) Schur canonical form, i.e. A, B are upper quasi
     !> triangular and D, E are upper triangular.
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve  Zx = scale b, where
     !> Z is defined as
     !> Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
     !> [ kron(In, D)  -kron(E**T, Im) ].
     !> Here Ik is the identity matrix of size k and X**T is the transpose of
     !> X. kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'T', QTGSYL solves the transposed system Z**T*y = scale*b,
     !> which is equivalent to solve for R and L in
     !> A**T * R + D**T * L = scale * C           (3)
     !> R * B**T + L * E**T = scale * -F
     !> This case (TRANS = 'T') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using QLACON.
     !> If IJOB >= 1, QTGSYL computes a Frobenius norm-based estimate
     !> of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z. See [1-2] for more
     !> information.
     !> This is a level 3 BLAS algorithm.

     pure subroutine la_qtgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(qp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           real(qp),intent(inout) :: c(ldc,*),f(ldf,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_qcopy by calls to la_qlaset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,ppqq,pq,q
           real(qp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: real,max,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine optimal block sizes mb and nb
           mb = la_ilaenv(2,'QTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'QTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_qlaset('F',m,n,zero,zero,c,ldc)
                 call la_qlaset('F',m,n,zero,zero,f,ldf)
              else if (ijob >= 1) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              loop_30: do iround = 1,isolve
                 ! use unblocked level 2 solver
                 dscale = zero
                 dsum = one
                 pq = 0
                 call la_qtgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,iwork,pq,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=qp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=qp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_qlacpy('F',m,n,c,ldc,work,m)
                    call la_qlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_qlaset('F',m,n,zero,zero,c,ldc)
                    call la_qlaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_qlacpy('F',m,n,work,m,c,ldc)
                    call la_qlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           if (a(i,i - 1) /= zero) i = i + 1
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           if (b(j,j - 1) /= zero) j = j + 1
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j)-subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1,..., 1; j = 1, 2,..., q
                 dscale = zero
                 dsum = one
                 pq = 0
                 scale = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       ppqq = 0
                       call la_qtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,iwork(q + 2),ppqq,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + ppqq
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_qscal(m,scaloc,c(1,k),1)
                             call la_qscal(m,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_qscal(is - 1,scaloc,c(1,k),1)
                             call la_qscal(is - 1,scaloc,f(1,k),1)
                          end do
                          do k = js,je
                             call la_qscal(m - ie,scaloc,c(ie + 1,k),1)
                             call la_qscal(m - ie,scaloc,f(ie + 1,k),1)
                          end do
                          do k = je + 1,n
                             call la_qscal(m,scaloc,c(1,k),1)
                             call la_qscal(m,scaloc,f(1,k),1)
                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i, j) and l(i, j) into remaining
                       ! equation.
                       if (i > 1) then
                          call la_qgemm('N','N',is - 1,nb,mb,-one,a(1,is),lda,c(is, &
                                    js),ldc,one,c(1,js),ldc)
                          call la_qgemm('N','N',is - 1,nb,mb,-one,d(1,is),ldd,c(is, &
                                    js),ldc,one,f(1,js),ldf)
                       end if
                       if (j < q) then
                          call la_qgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,b(js, &
                                    je + 1),ldb,one,c(is,je + 1),ldc)
                          call la_qgemm('N','N',mb,n - je,nb,one,f(is,js),ldf,e(js, &
                                    je + 1),lde,one,f(is,je + 1),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=qp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=qp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_qlacpy('F',m,n,c,ldc,work,m)
                    call la_qlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_qlaset('F',m,n,zero,zero,c,ldc)
                    call la_qlaset('F',m,n,zero,zero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_qlacpy('F',m,n,work,m,c,ldc)
                    call la_qlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                   ! a(i, i)**t * r(i, j)  + d(i, i)**t * l(i, j)  =  c(i, j)
                   ! r(i, j)  * b(j, j)**t + l(i, j)  * e(j, j)**t = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_qtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,iwork(q + 2),ppqq,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_qscal(m,scaloc,c(1,k),1)
                          call la_qscal(m,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_qscal(is - 1,scaloc,c(1,k),1)
                          call la_qscal(is - 1,scaloc,f(1,k),1)
                       end do
                       do k = js,je
                          call la_qscal(m - ie,scaloc,c(ie + 1,k),1)
                          call la_qscal(m - ie,scaloc,f(ie + 1,k),1)
                       end do
                       do k = je + 1,n
                          call la_qscal(m,scaloc,c(1,k),1)
                          call la_qscal(m,scaloc,f(1,k),1)
                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (j > p + 2) then
                       call la_qgemm('N','T',mb,js - 1,nb,one,c(is,js),ldc,b(1,js) &
                                 ,ldb,one,f(is,1),ldf)
                       call la_qgemm('N','T',mb,js - 1,nb,one,f(is,js),ldf,e(1,js) &
                                 ,lde,one,f(is,1),ldf)
                    end if
                    if (i < p) then
                       call la_qgemm('T','N',m - ie,nb,mb,-one,a(is,ie + 1),lda,c(is, &
                                 js),ldc,one,c(ie + 1,js),ldc)
                       call la_qgemm('T','N',m - ie,nb,mb,-one,d(is,ie + 1),ldd,f(is, &
                                 js),ldf,one,c(ie + 1,js),ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_qtgsyl
#endif

     !> SLAGV2: computes the Generalized Schur factorization of a real 2-by-2
     !> matrix pencil (A,B) where B is upper triangular. This routine
     !> computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
     !> SNR such that
     !> 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
     !> types), then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],
     !> 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
     !> then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]
     !> where b11 >= b22 > 0.

     pure subroutine la_slagv2(a,lda,b,ldb,alphar,alphai,beta,csl,snl,csr,snr)
        use la_constants_sp,only:zero,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb
           real(sp),intent(out) :: csl,csr,snl,snr
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(2),alphar(2),beta(2)
        ! =====================================================================

           ! Local Scalars
           real(sp) :: anorm,ascale,bnorm,bscale,h1,h2,h3,qq,r,rr,safmin,scale1, &
                     scale2,t,ulp,wi,wr1,wr2
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           safmin = la_slamch('S')
           ulp = la_slamch('P')
           ! scale a
           anorm = max(abs(a(1,1)) + abs(a(2,1)),abs(a(1,2)) + abs(a(2,2)), &
                     safmin)
           ascale = one/anorm
           a(1,1) = ascale*a(1,1)
           a(1,2) = ascale*a(1,2)
           a(2,1) = ascale*a(2,1)
           a(2,2) = ascale*a(2,2)
           ! scale b
           bnorm = max(abs(b(1,1)),abs(b(1,2)) + abs(b(2,2)),safmin)
           bscale = one/bnorm
           b(1,1) = bscale*b(1,1)
           b(1,2) = bscale*b(1,2)
           b(2,2) = bscale*b(2,2)
           ! check if a can be deflated
           if (abs(a(2,1)) <= ulp) then
              csl = one
              snl = zero
              csr = one
              snr = zero
              a(2,1) = zero
              b(2,1) = zero
              wi = zero
           ! check if b is singular
           else if (abs(b(1,1)) <= ulp) then
              call la_slartg(a(1,1),a(2,1),csl,snl,r)
              csr = one
              snr = zero
              call la_srot(2,a(1,1),lda,a(2,1),lda,csl,snl)
              call la_srot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
              a(2,1) = zero
              b(1,1) = zero
              b(2,1) = zero
              wi = zero
           else if (abs(b(2,2)) <= ulp) then
              call la_slartg(a(2,2),a(2,1),csr,snr,t)
              snr = -snr
              call la_srot(2,a(1,1),1,a(1,2),1,csr,snr)
              call la_srot(2,b(1,1),1,b(1,2),1,csr,snr)
              csl = one
              snl = zero
              a(2,1) = zero
              b(2,1) = zero
              b(2,2) = zero
              wi = zero
           else
              ! b is nonsingular, first compute the eigenvalues of (a,b)
              call la_slag2(a,lda,b,ldb,safmin,scale1,scale2,wr1,wr2,wi)
              if (wi == zero) then
                 ! two real eigenvalues, compute s*a-w*b
                 h1 = scale1*a(1,1) - wr1*b(1,1)
                 h2 = scale1*a(1,2) - wr1*b(1,2)
                 h3 = scale1*a(2,2) - wr1*b(2,2)
                 rr = la_slapy2(h1,h2)
                 qq = la_slapy2(scale1*a(2,1),h3)
                 if (rr > qq) then
                    ! find right rotation matrix to zero 1,1 element of
                    ! (sa - wb)
                    call la_slartg(h2,h1,csr,snr,t)
                 else
                    ! find right rotation matrix to zero 2,1 element of
                    ! (sa - wb)
                    call la_slartg(h3,scale1*a(2,1),csr,snr,t)
                 end if
                 snr = -snr
                 call la_srot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_srot(2,b(1,1),1,b(1,2),1,csr,snr)
                 ! compute inf norms of a and b
                 h1 = max(abs(a(1,1)) + abs(a(1,2)),abs(a(2,1)) + abs(a(2,2)))

                 h2 = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2)))

                 if ((scale1*h1) >= abs(wr1)*h2) then
                    ! find left rotation matrix q to zero out b(2,1)
                    call la_slartg(b(1,1),b(2,1),csl,snl,r)
                 else
                    ! find left rotation matrix q to zero out a(2,1)
                    call la_slartg(a(1,1),a(2,1),csl,snl,r)
                 end if
                 call la_srot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_srot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 a(2,1) = zero
                 b(2,1) = zero
              else
                 ! a pair of complex conjugate eigenvalues
                 ! first compute the svd of the matrix b
                 call la_slasv2(b(1,1),b(1,2),b(2,2),r,t,snr,csr,snl,csl)

                 ! form (a,b) := q(a,b)z**t where q is left rotation matrix and
                 ! z is right rotation matrix computed from la_slasv2
                 call la_srot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_srot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 call la_srot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_srot(2,b(1,1),1,b(1,2),1,csr,snr)
                 b(2,1) = zero
                 b(1,2) = zero
              end if
           end if
           ! unscaling
           a(1,1) = anorm*a(1,1)
           a(2,1) = anorm*a(2,1)
           a(1,2) = anorm*a(1,2)
           a(2,2) = anorm*a(2,2)
           b(1,1) = bnorm*b(1,1)
           b(2,1) = bnorm*b(2,1)
           b(1,2) = bnorm*b(1,2)
           b(2,2) = bnorm*b(2,2)
           if (wi == zero) then
              alphar(1) = a(1,1)
              alphar(2) = a(2,2)
              alphai(1) = zero
              alphai(2) = zero
              beta(1) = b(1,1)
              beta(2) = b(2,2)
           else
              alphar(1) = anorm*wr1/scale1/bnorm
              alphai(1) = anorm*wi/scale1/bnorm
              alphar(2) = alphar(1)
              alphai(2) = -alphai(1)
              beta(1) = one
              beta(2) = one
           end if
           return
     end subroutine la_slagv2
     !> DLAGV2: computes the Generalized Schur factorization of a real 2-by-2
     !> matrix pencil (A,B) where B is upper triangular. This routine
     !> computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
     !> SNR such that
     !> 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
     !> types), then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],
     !> 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
     !> then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]
     !> where b11 >= b22 > 0.

     pure subroutine la_dlagv2(a,lda,b,ldb,alphar,alphai,beta,csl,snl,csr,snr)
        use la_constants_dp,only:zero,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb
           real(dp),intent(out) :: csl,csr,snl,snr
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(2),alphar(2),beta(2)
        ! =====================================================================

           ! Local Scalars
           real(dp) :: anorm,ascale,bnorm,bscale,h1,h2,h3,qq,r,rr,safmin,scale1, &
                     scale2,t,ulp,wi,wr1,wr2
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           safmin = la_dlamch('S')
           ulp = la_dlamch('P')
           ! scale a
           anorm = max(abs(a(1,1)) + abs(a(2,1)),abs(a(1,2)) + abs(a(2,2)), &
                     safmin)
           ascale = one/anorm
           a(1,1) = ascale*a(1,1)
           a(1,2) = ascale*a(1,2)
           a(2,1) = ascale*a(2,1)
           a(2,2) = ascale*a(2,2)
           ! scale b
           bnorm = max(abs(b(1,1)),abs(b(1,2)) + abs(b(2,2)),safmin)
           bscale = one/bnorm
           b(1,1) = bscale*b(1,1)
           b(1,2) = bscale*b(1,2)
           b(2,2) = bscale*b(2,2)
           ! check if a can be deflated
           if (abs(a(2,1)) <= ulp) then
              csl = one
              snl = zero
              csr = one
              snr = zero
              a(2,1) = zero
              b(2,1) = zero
              wi = zero
           ! check if b is singular
           else if (abs(b(1,1)) <= ulp) then
              call la_dlartg(a(1,1),a(2,1),csl,snl,r)
              csr = one
              snr = zero
              call la_drot(2,a(1,1),lda,a(2,1),lda,csl,snl)
              call la_drot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
              a(2,1) = zero
              b(1,1) = zero
              b(2,1) = zero
              wi = zero
           else if (abs(b(2,2)) <= ulp) then
              call la_dlartg(a(2,2),a(2,1),csr,snr,t)
              snr = -snr
              call la_drot(2,a(1,1),1,a(1,2),1,csr,snr)
              call la_drot(2,b(1,1),1,b(1,2),1,csr,snr)
              csl = one
              snl = zero
              a(2,1) = zero
              b(2,1) = zero
              b(2,2) = zero
              wi = zero
           else
              ! b is nonsingular, first compute the eigenvalues of (a,b)
              call la_dlag2(a,lda,b,ldb,safmin,scale1,scale2,wr1,wr2,wi)
              if (wi == zero) then
                 ! two real eigenvalues, compute s*a-w*b
                 h1 = scale1*a(1,1) - wr1*b(1,1)
                 h2 = scale1*a(1,2) - wr1*b(1,2)
                 h3 = scale1*a(2,2) - wr1*b(2,2)
                 rr = la_dlapy2(h1,h2)
                 qq = la_dlapy2(scale1*a(2,1),h3)
                 if (rr > qq) then
                    ! find right rotation matrix to zero 1,1 element of
                    ! (sa - wb)
                    call la_dlartg(h2,h1,csr,snr,t)
                 else
                    ! find right rotation matrix to zero 2,1 element of
                    ! (sa - wb)
                    call la_dlartg(h3,scale1*a(2,1),csr,snr,t)
                 end if
                 snr = -snr
                 call la_drot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_drot(2,b(1,1),1,b(1,2),1,csr,snr)
                 ! compute inf norms of a and b
                 h1 = max(abs(a(1,1)) + abs(a(1,2)),abs(a(2,1)) + abs(a(2,2)))

                 h2 = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2)))

                 if ((scale1*h1) >= abs(wr1)*h2) then
                    ! find left rotation matrix q to zero out b(2,1)
                    call la_dlartg(b(1,1),b(2,1),csl,snl,r)
                 else
                    ! find left rotation matrix q to zero out a(2,1)
                    call la_dlartg(a(1,1),a(2,1),csl,snl,r)
                 end if
                 call la_drot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_drot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 a(2,1) = zero
                 b(2,1) = zero
              else
                 ! a pair of complex conjugate eigenvalues
                 ! first compute the svd of the matrix b
                 call la_dlasv2(b(1,1),b(1,2),b(2,2),r,t,snr,csr,snl,csl)

                 ! form (a,b) := q(a,b)z**t where q is left rotation matrix and
                 ! z is right rotation matrix computed from la_dlasv2
                 call la_drot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_drot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 call la_drot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_drot(2,b(1,1),1,b(1,2),1,csr,snr)
                 b(2,1) = zero
                 b(1,2) = zero
              end if
           end if
           ! unscaling
           a(1,1) = anorm*a(1,1)
           a(2,1) = anorm*a(2,1)
           a(1,2) = anorm*a(1,2)
           a(2,2) = anorm*a(2,2)
           b(1,1) = bnorm*b(1,1)
           b(2,1) = bnorm*b(2,1)
           b(1,2) = bnorm*b(1,2)
           b(2,2) = bnorm*b(2,2)
           if (wi == zero) then
              alphar(1) = a(1,1)
              alphar(2) = a(2,2)
              alphai(1) = zero
              alphai(2) = zero
              beta(1) = b(1,1)
              beta(2) = b(2,2)
           else
              alphar(1) = anorm*wr1/scale1/bnorm
              alphai(1) = anorm*wi/scale1/bnorm
              alphar(2) = alphar(1)
              alphai(2) = -alphai(1)
              beta(1) = one
              beta(2) = one
           end if
           return
     end subroutine la_dlagv2
#ifdef LA_WITH_XDP
     !> XLAGV2: computes the Generalized Schur factorization of a real 2-by-2
     !> matrix pencil (A,B) where B is upper triangular. This routine
     !> computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
     !> SNR such that
     !> 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
     !> types), then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],
     !> 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
     !> then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]
     !> where b11 >= b22 > 0.

     pure subroutine la_xlagv2(a,lda,b,ldb,alphar,alphai,beta,csl,snl,csr,snr)
        use la_constants_xdp,only:zero,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb
           real(xdp),intent(out) :: csl,csr,snl,snr
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: alphai(2),alphar(2),beta(2)
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: anorm,ascale,bnorm,bscale,h1,h2,h3,qq,r,rr,safmin,scale1, &
                     scale2,t,ulp,wi,wr1,wr2
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           safmin = la_xlamch('S')
           ulp = la_xlamch('P')
           ! scale a
           anorm = max(abs(a(1,1)) + abs(a(2,1)),abs(a(1,2)) + abs(a(2,2)), &
                     safmin)
           ascale = one/anorm
           a(1,1) = ascale*a(1,1)
           a(1,2) = ascale*a(1,2)
           a(2,1) = ascale*a(2,1)
           a(2,2) = ascale*a(2,2)
           ! scale b
           bnorm = max(abs(b(1,1)),abs(b(1,2)) + abs(b(2,2)),safmin)
           bscale = one/bnorm
           b(1,1) = bscale*b(1,1)
           b(1,2) = bscale*b(1,2)
           b(2,2) = bscale*b(2,2)
           ! check if a can be deflated
           if (abs(a(2,1)) <= ulp) then
              csl = one
              snl = zero
              csr = one
              snr = zero
              a(2,1) = zero
              b(2,1) = zero
              wi = zero
           ! check if b is singular
           else if (abs(b(1,1)) <= ulp) then
              call la_xlartg(a(1,1),a(2,1),csl,snl,r)
              csr = one
              snr = zero
              call la_xrot(2,a(1,1),lda,a(2,1),lda,csl,snl)
              call la_xrot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
              a(2,1) = zero
              b(1,1) = zero
              b(2,1) = zero
              wi = zero
           else if (abs(b(2,2)) <= ulp) then
              call la_xlartg(a(2,2),a(2,1),csr,snr,t)
              snr = -snr
              call la_xrot(2,a(1,1),1,a(1,2),1,csr,snr)
              call la_xrot(2,b(1,1),1,b(1,2),1,csr,snr)
              csl = one
              snl = zero
              a(2,1) = zero
              b(2,1) = zero
              b(2,2) = zero
              wi = zero
           else
              ! b is nonsingular, first compute the eigenvalues of (a,b)
              call la_xlag2(a,lda,b,ldb,safmin,scale1,scale2,wr1,wr2,wi)
              if (wi == zero) then
                 ! two real eigenvalues, compute s*a-w*b
                 h1 = scale1*a(1,1) - wr1*b(1,1)
                 h2 = scale1*a(1,2) - wr1*b(1,2)
                 h3 = scale1*a(2,2) - wr1*b(2,2)
                 rr = la_xlapy2(h1,h2)
                 qq = la_xlapy2(scale1*a(2,1),h3)
                 if (rr > qq) then
                    ! find right rotation matrix to zero 1,1 element of
                    ! (sa - wb)
                    call la_xlartg(h2,h1,csr,snr,t)
                 else
                    ! find right rotation matrix to zero 2,1 element of
                    ! (sa - wb)
                    call la_xlartg(h3,scale1*a(2,1),csr,snr,t)
                 end if
                 snr = -snr
                 call la_xrot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_xrot(2,b(1,1),1,b(1,2),1,csr,snr)
                 ! compute inf norms of a and b
                 h1 = max(abs(a(1,1)) + abs(a(1,2)),abs(a(2,1)) + abs(a(2,2)))

                 h2 = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2)))

                 if ((scale1*h1) >= abs(wr1)*h2) then
                    ! find left rotation matrix q to zero out b(2,1)
                    call la_xlartg(b(1,1),b(2,1),csl,snl,r)
                 else
                    ! find left rotation matrix q to zero out a(2,1)
                    call la_xlartg(a(1,1),a(2,1),csl,snl,r)
                 end if
                 call la_xrot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_xrot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 a(2,1) = zero
                 b(2,1) = zero
              else
                 ! a pair of complex conjugate eigenvalues
                 ! first compute the svd of the matrix b
                 call la_xlasv2(b(1,1),b(1,2),b(2,2),r,t,snr,csr,snl,csl)

                 ! form (a,b) := q(a,b)z**t where q is left rotation matrix and
                 ! z is right rotation matrix computed from la_xlasv2
                 call la_xrot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_xrot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 call la_xrot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_xrot(2,b(1,1),1,b(1,2),1,csr,snr)
                 b(2,1) = zero
                 b(1,2) = zero
              end if
           end if
           ! unscaling
           a(1,1) = anorm*a(1,1)
           a(2,1) = anorm*a(2,1)
           a(1,2) = anorm*a(1,2)
           a(2,2) = anorm*a(2,2)
           b(1,1) = bnorm*b(1,1)
           b(2,1) = bnorm*b(2,1)
           b(1,2) = bnorm*b(1,2)
           b(2,2) = bnorm*b(2,2)
           if (wi == zero) then
              alphar(1) = a(1,1)
              alphar(2) = a(2,2)
              alphai(1) = zero
              alphai(2) = zero
              beta(1) = b(1,1)
              beta(2) = b(2,2)
           else
              alphar(1) = anorm*wr1/scale1/bnorm
              alphai(1) = anorm*wi/scale1/bnorm
              alphar(2) = alphar(1)
              alphai(2) = -alphai(1)
              beta(1) = one
              beta(2) = one
           end if
           return
     end subroutine la_xlagv2
#endif
#ifdef LA_WITH_QP
     !> QLAGV2: computes the Generalized Schur factorization of a real 2-by-2
     !> matrix pencil (A,B) where B is upper triangular. This routine
     !> computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
     !> SNR such that
     !> 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
     !> types), then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],
     !> 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
     !> then
     !> [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
     !> [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
     !> [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
     !> [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]
     !> where b11 >= b22 > 0.

     pure subroutine la_qlagv2(a,lda,b,ldb,alphar,alphai,beta,csl,snl,csr,snr)
        use la_constants_qp,only:zero,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb
           real(qp),intent(out) :: csl,csr,snl,snr
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(2),alphar(2),beta(2)
        ! =====================================================================

           ! Local Scalars
           real(qp) :: anorm,ascale,bnorm,bscale,h1,h2,h3,qq,r,rr,safmin,scale1, &
                     scale2,t,ulp,wi,wr1,wr2
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           safmin = la_qlamch('S')
           ulp = la_qlamch('P')
           ! scale a
           anorm = max(abs(a(1,1)) + abs(a(2,1)),abs(a(1,2)) + abs(a(2,2)), &
                     safmin)
           ascale = one/anorm
           a(1,1) = ascale*a(1,1)
           a(1,2) = ascale*a(1,2)
           a(2,1) = ascale*a(2,1)
           a(2,2) = ascale*a(2,2)
           ! scale b
           bnorm = max(abs(b(1,1)),abs(b(1,2)) + abs(b(2,2)),safmin)
           bscale = one/bnorm
           b(1,1) = bscale*b(1,1)
           b(1,2) = bscale*b(1,2)
           b(2,2) = bscale*b(2,2)
           ! check if a can be deflated
           if (abs(a(2,1)) <= ulp) then
              csl = one
              snl = zero
              csr = one
              snr = zero
              a(2,1) = zero
              b(2,1) = zero
              wi = zero
           ! check if b is singular
           else if (abs(b(1,1)) <= ulp) then
              call la_qlartg(a(1,1),a(2,1),csl,snl,r)
              csr = one
              snr = zero
              call la_qrot(2,a(1,1),lda,a(2,1),lda,csl,snl)
              call la_qrot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
              a(2,1) = zero
              b(1,1) = zero
              b(2,1) = zero
              wi = zero
           else if (abs(b(2,2)) <= ulp) then
              call la_qlartg(a(2,2),a(2,1),csr,snr,t)
              snr = -snr
              call la_qrot(2,a(1,1),1,a(1,2),1,csr,snr)
              call la_qrot(2,b(1,1),1,b(1,2),1,csr,snr)
              csl = one
              snl = zero
              a(2,1) = zero
              b(2,1) = zero
              b(2,2) = zero
              wi = zero
           else
              ! b is nonsingular, first compute the eigenvalues of (a,b)
              call la_qlag2(a,lda,b,ldb,safmin,scale1,scale2,wr1,wr2,wi)
              if (wi == zero) then
                 ! two real eigenvalues, compute s*a-w*b
                 h1 = scale1*a(1,1) - wr1*b(1,1)
                 h2 = scale1*a(1,2) - wr1*b(1,2)
                 h3 = scale1*a(2,2) - wr1*b(2,2)
                 rr = la_qlapy2(h1,h2)
                 qq = la_qlapy2(scale1*a(2,1),h3)
                 if (rr > qq) then
                    ! find right rotation matrix to zero 1,1 element of
                    ! (sa - wb)
                    call la_qlartg(h2,h1,csr,snr,t)
                 else
                    ! find right rotation matrix to zero 2,1 element of
                    ! (sa - wb)
                    call la_qlartg(h3,scale1*a(2,1),csr,snr,t)
                 end if
                 snr = -snr
                 call la_qrot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_qrot(2,b(1,1),1,b(1,2),1,csr,snr)
                 ! compute inf norms of a and b
                 h1 = max(abs(a(1,1)) + abs(a(1,2)),abs(a(2,1)) + abs(a(2,2)))

                 h2 = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2)))

                 if ((scale1*h1) >= abs(wr1)*h2) then
                    ! find left rotation matrix q to zero out b(2,1)
                    call la_qlartg(b(1,1),b(2,1),csl,snl,r)
                 else
                    ! find left rotation matrix q to zero out a(2,1)
                    call la_qlartg(a(1,1),a(2,1),csl,snl,r)
                 end if
                 call la_qrot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_qrot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 a(2,1) = zero
                 b(2,1) = zero
              else
                 ! a pair of complex conjugate eigenvalues
                 ! first compute the svd of the matrix b
                 call la_qlasv2(b(1,1),b(1,2),b(2,2),r,t,snr,csr,snl,csl)

                 ! form (a,b) := q(a,b)z**t where q is left rotation matrix and
                 ! z is right rotation matrix computed from la_qlasv2
                 call la_qrot(2,a(1,1),lda,a(2,1),lda,csl,snl)
                 call la_qrot(2,b(1,1),ldb,b(2,1),ldb,csl,snl)
                 call la_qrot(2,a(1,1),1,a(1,2),1,csr,snr)
                 call la_qrot(2,b(1,1),1,b(1,2),1,csr,snr)
                 b(2,1) = zero
                 b(1,2) = zero
              end if
           end if
           ! unscaling
           a(1,1) = anorm*a(1,1)
           a(2,1) = anorm*a(2,1)
           a(1,2) = anorm*a(1,2)
           a(2,2) = anorm*a(2,2)
           b(1,1) = bnorm*b(1,1)
           b(2,1) = bnorm*b(2,1)
           b(1,2) = bnorm*b(1,2)
           b(2,2) = bnorm*b(2,2)
           if (wi == zero) then
              alphar(1) = a(1,1)
              alphar(2) = a(2,2)
              alphai(1) = zero
              alphai(2) = zero
              beta(1) = b(1,1)
              beta(2) = b(2,2)
           else
              alphar(1) = anorm*wr1/scale1/bnorm
              alphai(1) = anorm*wi/scale1/bnorm
              alphar(2) = alphar(1)
              alphai(2) = -alphai(1)
              beta(1) = one
              beta(2) = one
           end if
           return
     end subroutine la_qlagv2
#endif

     !> STGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of real matrices (S,P), where S is a quasi-triangular matrix
     !> and P is upper triangular.  Matrix pairs of this type are produced by
     !> the generalized Schur factorization of a matrix pair (A,B):
     !> A = Q*S*Z**T,  B = Q*P*Z**T
     !> as computed by SGGHRD + SHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal blocks of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the orthogonal factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_stgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(sp),intent(in) :: p(ldp,*),s(lds,*)
           real(sp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: safety = 1.0e+2_sp

           ! Local Scalars
           logical(lk) :: compl,compr,il2by2,ilabad,ilall,ilback,ilbbad,ilcomp,ilcplx, &
                     lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,iinfo,im,iside,j,ja,jc,je,jr,jw, &
                     na,nw
           real(sp) :: acoef,acoefa,anorm,ascale,bcoefa,bcoefi,bcoefr,big,bignum,bnorm, &
           bscale,cim2a,cim2b,cimaga,cimagb,cre2a,cre2b,creala,crealb,dmin,safmin, &
                     salfar,sbeta,scale,small,temp,temp2,temp2i,temp2r,ulp,xmax,xscale
           ! Local Arrays
           real(sp) :: bdiag(2),sum(2,2),sums(2,2),sump(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
              ilall = .true.
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('STGEVC',-info)
              return
           end if
           ! count the number of eigenvectors to be computed
           if (.not. ilall) then
              im = 0
              ilcplx = .false.
              loop_10: do j = 1,n
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_10
                 end if
                 if (j < n) then
                    if (s(j + 1,j) /= zero) ilcplx = .true.
                 end if
                 if (ilcplx) then
                    if (select(j) .or. select(j + 1)) im = im + 2
                 else
                    if (select(j)) im = im + 1
                 end if
              end do loop_10
           else
              im = n
           end if
           ! check 2-by-2 diagonal blocks of a, b
           ilabad = .false.
           ilbbad = .false.
           do j = 1,n - 1
              if (s(j + 1,j) /= zero) then
                 if (p(j,j) == zero .or. p(j + 1,j + 1) == zero .or. p(j,j + 1) /= zero) ilbbad = &
                           .true.
                 if (j < n - 1) then
                    if (s(j + 2,j + 1) /= zero) ilabad = .true.
                 end if
              end if
           end do
           if (ilabad) then
              info = -5
           else if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('STGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_slamch('SAFE MINIMUM')
           big = one/safmin
           call la_slabad(safmin,big)
           ulp = la_slamch('EPSILON')*la_slamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part (i.e., excluding all elements belonging to the diagonal
           ! blocks) of a and b to check for possible overflow in the
           ! triangular solver.
           anorm = abs(s(1,1))
           if (n > 1) anorm = anorm + abs(s(2,1))
           bnorm = abs(p(1,1))
           work(1) = zero
           work(n + 1) = zero
           do j = 2,n
              temp = zero
              temp2 = zero
              if (s(j,j - 1) == zero) then
                 iend = j - 1
              else
                 iend = j - 2
              end if
              do i = 1,iend
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              work(j) = temp
              work(n + j) = temp2
              do i = iend + 1,min(j + 1,n)
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              anorm = max(anorm,temp)
              bnorm = max(bnorm,temp2)
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_220: do je = 1,n
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at.
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_220
                 end if
                 nw = 1
                 if (je < n) then
                    if (s(je + 1,je) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je + 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_220
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       ieig = ieig + 1
                       do jr = 1,n
                          vl(jr,ieig) = zero
                       end do
                       vl(ieig,ieig) = one
                       cycle loop_220
                    end if
                 end if
                 ! clear vector
                 do jr = 1,nw*n
                    work(2*n + jr) = zero
                 end do
                                                       ! t
                 ! compute coefficients in  ( a a - b b )  y = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                 else
                    ! complex eigenvalue
                    call la_slag2(s(je,je),lds,p(je,je),ldp,safmin*safety,acoef, &
                              temp,bcoefr,temp2,bcoefi)
                    bcoefi = -bcoefi
                    if (bcoefi == zero) then
                       info = je
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    temp = acoef*s(je + 1,je)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) > abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je + 1) = -temp2r/temp
                       work(3*n + je + 1) = -temp2i/temp
                    else
                       work(2*n + je + 1) = one
                       work(3*n + je + 1) = zero
                       temp = acoef*s(je,je + 1)
                       work(2*n + je) = (bcoefr*p(je + 1,je + 1) - acoef*s(je + 1,je + 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je + 1,je + 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je + 1) &
                              ) + abs(work(3*n + je + 1)))
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                 ! t
                 ! triangular solve of  (a a - b b)  y = 0
                                         ! t
                 ! (rowwise in  (a a - b b) , or columnwise in (a a - b b) )
                 il2by2 = .false.
                 loop_160: do j = je + nw,n
                    if (il2by2) then
                       il2by2 = .false.
                       cycle loop_160
                    end if
                    na = 1
                    bdiag(1) = p(j,j)
                    if (j < n) then
                       if (s(j + 1,j) /= zero) then
                          il2by2 = .true.
                          bdiag(2) = p(j + 1,j + 1)
                          na = 2
                       end if
                    end if
                    ! check whether scaling is necessary for dot products
                    xscale = one/max(one,xmax)
                    temp = max(work(j),work(n + j),acoefa*work(j) + bcoefa*work(n + j))

                    if (il2by2) temp = max(temp,work(j + 1),work(n + j + 1),acoefa*work(j + 1) + &
                              bcoefa*work(n + j + 1))
                    if (temp > bignum*xscale) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = xmax*xscale
                    end if
                    ! compute dot products
                          ! j-1
                    ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                          ! k=je
                    ! to reduce the op count, this is done as
                    ! _        j-1                  _        j-1
                    ! a*conjg( sum  s(k,j)*x(k) ) - b*conjg( sum  p(k,j)*x(k) )
                             ! k=je                          k=je
                    ! which may cause underflow problems if a or b are close
                    ! to underflow.  (e.g., less than small.)
                    do jw = 1,nw
                       do ja = 1,na
                          sums(ja,jw) = zero
                          sump(ja,jw) = zero
                          do jr = je,j - 1
                             sums(ja,jw) = sums(ja,jw) + s(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                             sump(ja,jw) = sump(ja,jw) + p(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                          end do
                       end do
                    end do
                    do ja = 1,na
                       if (ilcplx) then
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1) - bcoefi*sump( &
                                    ja,2)
                          sum(ja,2) = -acoef*sums(ja,2) + bcoefr*sump(ja,2) + bcoefi*sump( &
                                    ja,1)
                       else
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1)
                       end if
                    end do
                                        ! t
                    ! solve  ( a a - b b )  y = sum(,)
                    ! with scaling and perturbation of the denominator
                    call la_slaln2(.true.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),sum,2,bcoefr,bcoefi,work(2*n + j),n,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = scale*xmax
                    end if
                    xmax = max(xmax,temp)
                 end do loop_160
                 ! copy eigenvector to vl, back transforming if
                 ! howmny='b'.
                 ieig = ieig + 1
                 if (ilback) then
                    do jw = 0,nw - 1
                       call la_sgemv('N',n,n + 1 - je,one,vl(1,je),ldvl,work((jw + 2)*n + &
                                 je),1,zero,work((jw + 4)*n + 1),1)
                    end do
                    call la_slacpy(' ',n,nw,work(4*n + 1),n,vl(1,je),ldvl)
                    ibeg = 1
                 else
                    call la_slacpy(' ',n,nw,work(2*n + 1),n,vl(1,ieig),ldvl)
                    ibeg = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)) + abs(vl(j,ieig + 1)))
                    end do
                 else
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = ibeg,n
                          vl(jr,ieig + jw) = xscale*vl(jr,ieig + jw)
                       end do
                    end do
                 end if
                 ieig = ieig + nw - 1
              end do loop_220
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_500: do je = n,1,-1
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at -- if complex, select(je)
                 ! or select(je-1).
                 ! if this is a complex pair, the 2-by-2 diagonal block
                 ! corresponding to the eigenvalue is in rows/columns je-1:je
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_500
                 end if
                 nw = 1
                 if (je > 1) then
                    if (s(je,je - 1) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je - 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_500
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- unit eigenvector
                       ieig = ieig - 1
                       do jr = 1,n
                          vr(jr,ieig) = zero
                       end do
                       vr(ieig,ieig) = one
                       cycle loop_500
                    end if
                 end if
                 ! clear vector
                 do jw = 0,nw - 1
                    do jr = 1,n
                       work((jw + 2)*n + jr) = zero
                    end do
                 end do
                 ! compute coefficients in  ( a a - b b ) x = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                    ! compute contribution from column je of a and b to sum
                    ! (see "further details", above.)
                    do jr = 1,je - 1
                       work(2*n + jr) = bcoefr*p(jr,je) - acoef*s(jr,je)
                    end do
                 else
                    ! complex eigenvalue
                    call la_slag2(s(je - 1,je - 1),lds,p(je - 1,je - 1),ldp,safmin*safety, &
                              acoef,temp,bcoefr,temp2,bcoefi)
                    if (bcoefi == zero) then
                       info = je - 1
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    ! and contribution to sums
                    temp = acoef*s(je,je - 1)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) >= abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je - 1) = -temp2r/temp
                       work(3*n + je - 1) = -temp2i/temp
                    else
                       work(2*n + je - 1) = one
                       work(3*n + je - 1) = zero
                       temp = acoef*s(je - 1,je)
                       work(2*n + je) = (bcoefr*p(je - 1,je - 1) - acoef*s(je - 1,je - 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je - 1,je - 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je - 1) &
                              ) + abs(work(3*n + je - 1)))
                    ! compute contribution from columns je and je-1
                    ! of a and b to the sums.
                    creala = acoef*work(2*n + je - 1)
                    cimaga = acoef*work(3*n + je - 1)
                    crealb = bcoefr*work(2*n + je - 1) - bcoefi*work(3*n + je - 1)
                    cimagb = bcoefi*work(2*n + je - 1) + bcoefr*work(3*n + je - 1)
                    cre2a = acoef*work(2*n + je)
                    cim2a = acoef*work(3*n + je)
                    cre2b = bcoefr*work(2*n + je) - bcoefi*work(3*n + je)
                    cim2b = bcoefi*work(2*n + je) + bcoefr*work(3*n + je)
                    do jr = 1,je - 2
                       work(2*n + jr) = -creala*s(jr,je - 1) + crealb*p(jr,je - 1) - cre2a*s(jr, &
                                 je) + cre2b*p(jr,je)
                       work(3*n + jr) = -cimaga*s(jr,je - 1) + cimagb*p(jr,je - 1) - cim2a*s(jr, &
                                 je) + cim2b*p(jr,je)
                    end do
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                 ! columnwise triangular solve of  (a a - b b)  x = 0
                 il2by2 = .false.
                 loop_370: do j = je - nw,1,-1
                    ! if a 2-by-2 block, is in position j-1:j, wait until
                    ! next iteration to process it (when it will be j:j+1)
                    if (.not. il2by2 .and. j > 1) then
                       if (s(j,j - 1) /= zero) then
                          il2by2 = .true.
                          cycle loop_370
                       end if
                    end if
                    bdiag(1) = p(j,j)
                    if (il2by2) then
                       na = 2
                       bdiag(2) = p(j + 1,j + 1)
                    else
                       na = 1
                    end if
                    ! compute x(j) (and x(j+1), if 2-by-2 block)
                    call la_slaln2(.false.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),work(2*n + j),n,bcoefr,bcoefi,sum,2,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = 1,je
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                    end if
                    xmax = max(scale*xmax,temp)
                    do jw = 1,nw
                       do ja = 1,na
                          work((jw + 1)*n + j + ja - 1) = sum(ja,jw)
                       end do
                    end do
                    ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                    if (j > 1) then
                       ! check whether scaling is necessary for sum.
                       xscale = one/max(one,xmax)
                       temp = acoefa*work(j) + bcoefa*work(n + j)
                       if (il2by2) temp = max(temp,acoefa*work(j + 1) + bcoefa*work(n + j + 1))

                       temp = max(temp,acoefa,bcoefa)
                       if (temp > bignum*xscale) then
                          do jw = 0,nw - 1
                             do jr = 1,je
                                work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                             end do
                          end do
                          xmax = xmax*xscale
                       end if
                       ! compute the contributions of the off-diagonals of
                       ! column j (and j+1, if 2-by-2 block) of a and b to the
                       ! sums.
                       do ja = 1,na
                          if (ilcplx) then
                             creala = acoef*work(2*n + j + ja - 1)
                             cimaga = acoef*work(3*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1) - bcoefi*work(3*n + j + ja - 1)
                             cimagb = bcoefi*work(2*n + j + ja - 1) + bcoefr*work(3*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                                work(3*n + jr) = work(3*n + jr) - cimaga*s(jr,j + ja - 1) + cimagb*p( &
                                           jr,j + ja - 1)
                             end do
                          else
                             creala = acoef*work(2*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                             end do
                          end if
                       end do
                    end if
                    il2by2 = .false.
                 end do loop_370
                 ! copy eigenvector to vr, back transforming if
                 ! howmny='b'.
                 ieig = ieig - nw
                 if (ilback) then
                    do jw = 0,nw - 1
                       do jr = 1,n
                          work((jw + 4)*n + jr) = work((jw + 2)*n + 1)*vr(jr,1)
                       end do
                       ! a series of compiler directives to defeat
                       ! vectorization for the next loop
                       do jc = 2,je
                          do jr = 1,n
                             work((jw + 4)*n + jr) = work((jw + 4)*n + jr) + work((jw + 2)*n + jc) &
                                       *vr(jr,jc)
                          end do
                       end do
                    end do
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 4)*n + jr)
                       end do
                    end do
                    iend = n
                 else
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 2)*n + jr)
                       end do
                    end do
                    iend = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)) + abs(vr(j,ieig + 1)))
                    end do
                 else
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = 1,iend
                          vr(jr,ieig + jw) = xscale*vr(jr,ieig + jw)
                       end do
                    end do
                 end if
              end do loop_500
           end if
           return
     end subroutine la_stgevc
     !> DTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of real matrices (S,P), where S is a quasi-triangular matrix
     !> and P is upper triangular.  Matrix pairs of this type are produced by
     !> the generalized Schur factorization of a matrix pair (A,B):
     !> A = Q*S*Z**T,  B = Q*P*Z**T
     !> as computed by DGGHRD + DHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal blocks of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the orthogonal factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_dtgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(dp),intent(in) :: p(ldp,*),s(lds,*)
           real(dp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: safety = 1.0e+2_dp

           ! Local Scalars
           logical(lk) :: compl,compr,il2by2,ilabad,ilall,ilback,ilbbad,ilcomp,ilcplx, &
                     lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,iinfo,im,iside,j,ja,jc,je,jr,jw, &
                     na,nw
           real(dp) :: acoef,acoefa,anorm,ascale,bcoefa,bcoefi,bcoefr,big,bignum,bnorm, &
           bscale,cim2a,cim2b,cimaga,cimagb,cre2a,cre2b,creala,crealb,dmin,safmin, &
                     salfar,sbeta,scale,small,temp,temp2,temp2i,temp2r,ulp,xmax,xscale
           ! Local Arrays
           real(dp) :: bdiag(2),sum(2,2),sums(2,2),sump(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
              ilall = .true.
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors to be computed
           if (.not. ilall) then
              im = 0
              ilcplx = .false.
              loop_10: do j = 1,n
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_10
                 end if
                 if (j < n) then
                    if (s(j + 1,j) /= zero) ilcplx = .true.
                 end if
                 if (ilcplx) then
                    if (select(j) .or. select(j + 1)) im = im + 2
                 else
                    if (select(j)) im = im + 1
                 end if
              end do loop_10
           else
              im = n
           end if
           ! check 2-by-2 diagonal blocks of a, b
           ilabad = .false.
           ilbbad = .false.
           do j = 1,n - 1
              if (s(j + 1,j) /= zero) then
                 if (p(j,j) == zero .or. p(j + 1,j + 1) == zero .or. p(j,j + 1) /= zero) ilbbad = &
                           .true.
                 if (j < n - 1) then
                    if (s(j + 2,j + 1) /= zero) ilabad = .true.
                 end if
              end if
           end do
           if (ilabad) then
              info = -5
           else if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('DTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_dlamch('SAFE MINIMUM')
           big = one/safmin
           call la_dlabad(safmin,big)
           ulp = la_dlamch('EPSILON')*la_dlamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part (i.e., excluding all elements belonging to the diagonal
           ! blocks) of a and b to check for possible overflow in the
           ! triangular solver.
           anorm = abs(s(1,1))
           if (n > 1) anorm = anorm + abs(s(2,1))
           bnorm = abs(p(1,1))
           work(1) = zero
           work(n + 1) = zero
           do j = 2,n
              temp = zero
              temp2 = zero
              if (s(j,j - 1) == zero) then
                 iend = j - 1
              else
                 iend = j - 2
              end if
              do i = 1,iend
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              work(j) = temp
              work(n + j) = temp2
              do i = iend + 1,min(j + 1,n)
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              anorm = max(anorm,temp)
              bnorm = max(bnorm,temp2)
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_220: do je = 1,n
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at.
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_220
                 end if
                 nw = 1
                 if (je < n) then
                    if (s(je + 1,je) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je + 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_220
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       ieig = ieig + 1
                       do jr = 1,n
                          vl(jr,ieig) = zero
                       end do
                       vl(ieig,ieig) = one
                       cycle loop_220
                    end if
                 end if
                 ! clear vector
                 do jr = 1,nw*n
                    work(2*n + jr) = zero
                 end do
                                                       ! t
                 ! compute coefficients in  ( a a - b b )  y = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                 else
                    ! complex eigenvalue
                    call la_dlag2(s(je,je),lds,p(je,je),ldp,safmin*safety,acoef, &
                              temp,bcoefr,temp2,bcoefi)
                    bcoefi = -bcoefi
                    if (bcoefi == zero) then
                       info = je
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    temp = acoef*s(je + 1,je)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) > abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je + 1) = -temp2r/temp
                       work(3*n + je + 1) = -temp2i/temp
                    else
                       work(2*n + je + 1) = one
                       work(3*n + je + 1) = zero
                       temp = acoef*s(je,je + 1)
                       work(2*n + je) = (bcoefr*p(je + 1,je + 1) - acoef*s(je + 1,je + 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je + 1,je + 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je + 1) &
                              ) + abs(work(3*n + je + 1)))
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                 ! t
                 ! triangular solve of  (a a - b b)  y = 0
                                         ! t
                 ! (rowwise in  (a a - b b) , or columnwise in (a a - b b) )
                 il2by2 = .false.
                 loop_160: do j = je + nw,n
                    if (il2by2) then
                       il2by2 = .false.
                       cycle loop_160
                    end if
                    na = 1
                    bdiag(1) = p(j,j)
                    if (j < n) then
                       if (s(j + 1,j) /= zero) then
                          il2by2 = .true.
                          bdiag(2) = p(j + 1,j + 1)
                          na = 2
                       end if
                    end if
                    ! check whether scaling is necessary for dot products
                    xscale = one/max(one,xmax)
                    temp = max(work(j),work(n + j),acoefa*work(j) + bcoefa*work(n + j))

                    if (il2by2) temp = max(temp,work(j + 1),work(n + j + 1),acoefa*work(j + 1) + &
                              bcoefa*work(n + j + 1))
                    if (temp > bignum*xscale) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = xmax*xscale
                    end if
                    ! compute dot products
                          ! j-1
                    ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                          ! k=je
                    ! to reduce the op count, this is done as
                    ! _        j-1                  _        j-1
                    ! a*conjg( sum  s(k,j)*x(k) ) - b*conjg( sum  p(k,j)*x(k) )
                             ! k=je                          k=je
                    ! which may cause underflow problems if a or b are close
                    ! to underflow.  (e.g., less than small.)
                    do jw = 1,nw
                       do ja = 1,na
                          sums(ja,jw) = zero
                          sump(ja,jw) = zero
                          do jr = je,j - 1
                             sums(ja,jw) = sums(ja,jw) + s(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                             sump(ja,jw) = sump(ja,jw) + p(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                          end do
                       end do
                    end do
                    do ja = 1,na
                       if (ilcplx) then
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1) - bcoefi*sump( &
                                    ja,2)
                          sum(ja,2) = -acoef*sums(ja,2) + bcoefr*sump(ja,2) + bcoefi*sump( &
                                    ja,1)
                       else
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1)
                       end if
                    end do
                                        ! t
                    ! solve  ( a a - b b )  y = sum(,)
                    ! with scaling and perturbation of the denominator
                    call la_dlaln2(.true.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),sum,2,bcoefr,bcoefi,work(2*n + j),n,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = scale*xmax
                    end if
                    xmax = max(xmax,temp)
                 end do loop_160
                 ! copy eigenvector to vl, back transforming if
                 ! howmny='b'.
                 ieig = ieig + 1
                 if (ilback) then
                    do jw = 0,nw - 1
                       call la_dgemv('N',n,n + 1 - je,one,vl(1,je),ldvl,work((jw + 2)*n + &
                                 je),1,zero,work((jw + 4)*n + 1),1)
                    end do
                    call la_dlacpy(' ',n,nw,work(4*n + 1),n,vl(1,je),ldvl)
                    ibeg = 1
                 else
                    call la_dlacpy(' ',n,nw,work(2*n + 1),n,vl(1,ieig),ldvl)
                    ibeg = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)) + abs(vl(j,ieig + 1)))
                    end do
                 else
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = ibeg,n
                          vl(jr,ieig + jw) = xscale*vl(jr,ieig + jw)
                       end do
                    end do
                 end if
                 ieig = ieig + nw - 1
              end do loop_220
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_500: do je = n,1,-1
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at -- if complex, select(je)
                 ! or select(je-1).
                 ! if this is a complex pair, the 2-by-2 diagonal block
                 ! corresponding to the eigenvalue is in rows/columns je-1:je
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_500
                 end if
                 nw = 1
                 if (je > 1) then
                    if (s(je,je - 1) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je - 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_500
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- unit eigenvector
                       ieig = ieig - 1
                       do jr = 1,n
                          vr(jr,ieig) = zero
                       end do
                       vr(ieig,ieig) = one
                       cycle loop_500
                    end if
                 end if
                 ! clear vector
                 do jw = 0,nw - 1
                    do jr = 1,n
                       work((jw + 2)*n + jr) = zero
                    end do
                 end do
                 ! compute coefficients in  ( a a - b b ) x = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                    ! compute contribution from column je of a and b to sum
                    ! (see "further details", above.)
                    do jr = 1,je - 1
                       work(2*n + jr) = bcoefr*p(jr,je) - acoef*s(jr,je)
                    end do
                 else
                    ! complex eigenvalue
                    call la_dlag2(s(je - 1,je - 1),lds,p(je - 1,je - 1),ldp,safmin*safety, &
                              acoef,temp,bcoefr,temp2,bcoefi)
                    if (bcoefi == zero) then
                       info = je - 1
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    ! and contribution to sums
                    temp = acoef*s(je,je - 1)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) >= abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je - 1) = -temp2r/temp
                       work(3*n + je - 1) = -temp2i/temp
                    else
                       work(2*n + je - 1) = one
                       work(3*n + je - 1) = zero
                       temp = acoef*s(je - 1,je)
                       work(2*n + je) = (bcoefr*p(je - 1,je - 1) - acoef*s(je - 1,je - 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je - 1,je - 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je - 1) &
                              ) + abs(work(3*n + je - 1)))
                    ! compute contribution from columns je and je-1
                    ! of a and b to the sums.
                    creala = acoef*work(2*n + je - 1)
                    cimaga = acoef*work(3*n + je - 1)
                    crealb = bcoefr*work(2*n + je - 1) - bcoefi*work(3*n + je - 1)
                    cimagb = bcoefi*work(2*n + je - 1) + bcoefr*work(3*n + je - 1)
                    cre2a = acoef*work(2*n + je)
                    cim2a = acoef*work(3*n + je)
                    cre2b = bcoefr*work(2*n + je) - bcoefi*work(3*n + je)
                    cim2b = bcoefi*work(2*n + je) + bcoefr*work(3*n + je)
                    do jr = 1,je - 2
                       work(2*n + jr) = -creala*s(jr,je - 1) + crealb*p(jr,je - 1) - cre2a*s(jr, &
                                 je) + cre2b*p(jr,je)
                       work(3*n + jr) = -cimaga*s(jr,je - 1) + cimagb*p(jr,je - 1) - cim2a*s(jr, &
                                 je) + cim2b*p(jr,je)
                    end do
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                 ! columnwise triangular solve of  (a a - b b)  x = 0
                 il2by2 = .false.
                 loop_370: do j = je - nw,1,-1
                    ! if a 2-by-2 block, is in position j-1:j, wait until
                    ! next iteration to process it (when it will be j:j+1)
                    if (.not. il2by2 .and. j > 1) then
                       if (s(j,j - 1) /= zero) then
                          il2by2 = .true.
                          cycle loop_370
                       end if
                    end if
                    bdiag(1) = p(j,j)
                    if (il2by2) then
                       na = 2
                       bdiag(2) = p(j + 1,j + 1)
                    else
                       na = 1
                    end if
                    ! compute x(j) (and x(j+1), if 2-by-2 block)
                    call la_dlaln2(.false.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),work(2*n + j),n,bcoefr,bcoefi,sum,2,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = 1,je
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                    end if
                    xmax = max(scale*xmax,temp)
                    do jw = 1,nw
                       do ja = 1,na
                          work((jw + 1)*n + j + ja - 1) = sum(ja,jw)
                       end do
                    end do
                    ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                    if (j > 1) then
                       ! check whether scaling is necessary for sum.
                       xscale = one/max(one,xmax)
                       temp = acoefa*work(j) + bcoefa*work(n + j)
                       if (il2by2) temp = max(temp,acoefa*work(j + 1) + bcoefa*work(n + j + 1))

                       temp = max(temp,acoefa,bcoefa)
                       if (temp > bignum*xscale) then
                          do jw = 0,nw - 1
                             do jr = 1,je
                                work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                             end do
                          end do
                          xmax = xmax*xscale
                       end if
                       ! compute the contributions of the off-diagonals of
                       ! column j (and j+1, if 2-by-2 block) of a and b to the
                       ! sums.
                       do ja = 1,na
                          if (ilcplx) then
                             creala = acoef*work(2*n + j + ja - 1)
                             cimaga = acoef*work(3*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1) - bcoefi*work(3*n + j + ja - 1)
                             cimagb = bcoefi*work(2*n + j + ja - 1) + bcoefr*work(3*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                                work(3*n + jr) = work(3*n + jr) - cimaga*s(jr,j + ja - 1) + cimagb*p( &
                                           jr,j + ja - 1)
                             end do
                          else
                             creala = acoef*work(2*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                             end do
                          end if
                       end do
                    end if
                    il2by2 = .false.
                 end do loop_370
                 ! copy eigenvector to vr, back transforming if
                 ! howmny='b'.
                 ieig = ieig - nw
                 if (ilback) then
                    do jw = 0,nw - 1
                       do jr = 1,n
                          work((jw + 4)*n + jr) = work((jw + 2)*n + 1)*vr(jr,1)
                       end do
                       ! a series of compiler directives to defeat
                       ! vectorization for the next loop
                       do jc = 2,je
                          do jr = 1,n
                             work((jw + 4)*n + jr) = work((jw + 4)*n + jr) + work((jw + 2)*n + jc) &
                                       *vr(jr,jc)
                          end do
                       end do
                    end do
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 4)*n + jr)
                       end do
                    end do
                    iend = n
                 else
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 2)*n + jr)
                       end do
                    end do
                    iend = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)) + abs(vr(j,ieig + 1)))
                    end do
                 else
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = 1,iend
                          vr(jr,ieig + jw) = xscale*vr(jr,ieig + jw)
                       end do
                    end do
                 end if
              end do loop_500
           end if
           return
     end subroutine la_dtgevc
#ifdef LA_WITH_XDP
     !> XTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of real matrices (S,P), where S is a quasi-triangular matrix
     !> and P is upper triangular.  Matrix pairs of this type are produced by
     !> the generalized Schur factorization of a matrix pair (A,B):
     !> A = Q*S*Z**T,  B = Q*P*Z**T
     !> as computed by XGGHRD + XHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal blocks of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the orthogonal factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_xtgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(xdp),intent(in) :: p(ldp,*),s(lds,*)
           real(xdp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: safety = 1.0e+2_xdp

           ! Local Scalars
           logical(lk) :: compl,compr,il2by2,ilabad,ilall,ilback,ilbbad,ilcomp,ilcplx, &
                     lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,iinfo,im,iside,j,ja,jc,je,jr,jw, &
                     na,nw
           real(xdp) :: acoef,acoefa,anorm,ascale,bcoefa,bcoefi,bcoefr,big,bignum,bnorm, &
           bscale,cim2a,cim2b,cimaga,cimagb,cre2a,cre2b,creala,crealb,dmin,safmin, &
                     salfar,sbeta,scale,small,temp,temp2,temp2i,temp2r,ulp,xmax,xscale
           ! Local Arrays
           real(xdp) :: bdiag(2),sum(2,2),sums(2,2),sump(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
              ilall = .true.
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('XTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors to be computed
           if (.not. ilall) then
              im = 0
              ilcplx = .false.
              loop_10: do j = 1,n
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_10
                 end if
                 if (j < n) then
                    if (s(j + 1,j) /= zero) ilcplx = .true.
                 end if
                 if (ilcplx) then
                    if (select(j) .or. select(j + 1)) im = im + 2
                 else
                    if (select(j)) im = im + 1
                 end if
              end do loop_10
           else
              im = n
           end if
           ! check 2-by-2 diagonal blocks of a, b
           ilabad = .false.
           ilbbad = .false.
           do j = 1,n - 1
              if (s(j + 1,j) /= zero) then
                 if (p(j,j) == zero .or. p(j + 1,j + 1) == zero .or. p(j,j + 1) /= zero) ilbbad = &
                           .true.
                 if (j < n - 1) then
                    if (s(j + 2,j + 1) /= zero) ilabad = .true.
                 end if
              end if
           end do
           if (ilabad) then
              info = -5
           else if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('XTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_xlamch('SAFE MINIMUM')
           big = one/safmin
           call la_xlabad(safmin,big)
           ulp = la_xlamch('EPSILON')*la_xlamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part (i.e., excluding all elements belonging to the diagonal
           ! blocks) of a and b to check for possible overflow in the
           ! triangular solver.
           anorm = abs(s(1,1))
           if (n > 1) anorm = anorm + abs(s(2,1))
           bnorm = abs(p(1,1))
           work(1) = zero
           work(n + 1) = zero
           do j = 2,n
              temp = zero
              temp2 = zero
              if (s(j,j - 1) == zero) then
                 iend = j - 1
              else
                 iend = j - 2
              end if
              do i = 1,iend
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              work(j) = temp
              work(n + j) = temp2
              do i = iend + 1,min(j + 1,n)
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              anorm = max(anorm,temp)
              bnorm = max(bnorm,temp2)
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_220: do je = 1,n
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at.
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_220
                 end if
                 nw = 1
                 if (je < n) then
                    if (s(je + 1,je) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je + 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_220
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       ieig = ieig + 1
                       do jr = 1,n
                          vl(jr,ieig) = zero
                       end do
                       vl(ieig,ieig) = one
                       cycle loop_220
                    end if
                 end if
                 ! clear vector
                 do jr = 1,nw*n
                    work(2*n + jr) = zero
                 end do
                                                       ! t
                 ! compute coefficients in  ( a a - b b )  y = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                 else
                    ! complex eigenvalue
                    call la_xlag2(s(je,je),lds,p(je,je),ldp,safmin*safety,acoef, &
                              temp,bcoefr,temp2,bcoefi)
                    bcoefi = -bcoefi
                    if (bcoefi == zero) then
                       info = je
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    temp = acoef*s(je + 1,je)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) > abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je + 1) = -temp2r/temp
                       work(3*n + je + 1) = -temp2i/temp
                    else
                       work(2*n + je + 1) = one
                       work(3*n + je + 1) = zero
                       temp = acoef*s(je,je + 1)
                       work(2*n + je) = (bcoefr*p(je + 1,je + 1) - acoef*s(je + 1,je + 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je + 1,je + 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je + 1) &
                              ) + abs(work(3*n + je + 1)))
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                 ! t
                 ! triangular solve of  (a a - b b)  y = 0
                                         ! t
                 ! (rowwise in  (a a - b b) , or columnwise in (a a - b b) )
                 il2by2 = .false.
                 loop_160: do j = je + nw,n
                    if (il2by2) then
                       il2by2 = .false.
                       cycle loop_160
                    end if
                    na = 1
                    bdiag(1) = p(j,j)
                    if (j < n) then
                       if (s(j + 1,j) /= zero) then
                          il2by2 = .true.
                          bdiag(2) = p(j + 1,j + 1)
                          na = 2
                       end if
                    end if
                    ! check whether scaling is necessary for dot products
                    xscale = one/max(one,xmax)
                    temp = max(work(j),work(n + j),acoefa*work(j) + bcoefa*work(n + j))

                    if (il2by2) temp = max(temp,work(j + 1),work(n + j + 1),acoefa*work(j + 1) + &
                              bcoefa*work(n + j + 1))
                    if (temp > bignum*xscale) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = xmax*xscale
                    end if
                    ! compute dot products
                          ! j-1
                    ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                          ! k=je
                    ! to reduce the op count, this is done as
                    ! _        j-1                  _        j-1
                    ! a*conjg( sum  s(k,j)*x(k) ) - b*conjg( sum  p(k,j)*x(k) )
                             ! k=je                          k=je
                    ! which may cause underflow problems if a or b are close
                    ! to underflow.  (e.g., less than small.)
                    do jw = 1,nw
                       do ja = 1,na
                          sums(ja,jw) = zero
                          sump(ja,jw) = zero
                          do jr = je,j - 1
                             sums(ja,jw) = sums(ja,jw) + s(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                             sump(ja,jw) = sump(ja,jw) + p(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                          end do
                       end do
                    end do
                    do ja = 1,na
                       if (ilcplx) then
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1) - bcoefi*sump( &
                                    ja,2)
                          sum(ja,2) = -acoef*sums(ja,2) + bcoefr*sump(ja,2) + bcoefi*sump( &
                                    ja,1)
                       else
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1)
                       end if
                    end do
                                        ! t
                    ! solve  ( a a - b b )  y = sum(,)
                    ! with scaling and perturbation of the denominator
                    call la_xlaln2(.true.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),sum,2,bcoefr,bcoefi,work(2*n + j),n,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = scale*xmax
                    end if
                    xmax = max(xmax,temp)
                 end do loop_160
                 ! copy eigenvector to vl, back transforming if
                 ! howmny='b'.
                 ieig = ieig + 1
                 if (ilback) then
                    do jw = 0,nw - 1
                       call la_xgemv('N',n,n + 1 - je,one,vl(1,je),ldvl,work((jw + 2)*n + &
                                 je),1,zero,work((jw + 4)*n + 1),1)
                    end do
                    call la_xlacpy(' ',n,nw,work(4*n + 1),n,vl(1,je),ldvl)
                    ibeg = 1
                 else
                    call la_xlacpy(' ',n,nw,work(2*n + 1),n,vl(1,ieig),ldvl)
                    ibeg = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)) + abs(vl(j,ieig + 1)))
                    end do
                 else
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = ibeg,n
                          vl(jr,ieig + jw) = xscale*vl(jr,ieig + jw)
                       end do
                    end do
                 end if
                 ieig = ieig + nw - 1
              end do loop_220
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_500: do je = n,1,-1
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at -- if complex, select(je)
                 ! or select(je-1).
                 ! if this is a complex pair, the 2-by-2 diagonal block
                 ! corresponding to the eigenvalue is in rows/columns je-1:je
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_500
                 end if
                 nw = 1
                 if (je > 1) then
                    if (s(je,je - 1) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je - 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_500
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- unit eigenvector
                       ieig = ieig - 1
                       do jr = 1,n
                          vr(jr,ieig) = zero
                       end do
                       vr(ieig,ieig) = one
                       cycle loop_500
                    end if
                 end if
                 ! clear vector
                 do jw = 0,nw - 1
                    do jr = 1,n
                       work((jw + 2)*n + jr) = zero
                    end do
                 end do
                 ! compute coefficients in  ( a a - b b ) x = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                    ! compute contribution from column je of a and b to sum
                    ! (see "further details", above.)
                    do jr = 1,je - 1
                       work(2*n + jr) = bcoefr*p(jr,je) - acoef*s(jr,je)
                    end do
                 else
                    ! complex eigenvalue
                    call la_xlag2(s(je - 1,je - 1),lds,p(je - 1,je - 1),ldp,safmin*safety, &
                              acoef,temp,bcoefr,temp2,bcoefi)
                    if (bcoefi == zero) then
                       info = je - 1
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    ! and contribution to sums
                    temp = acoef*s(je,je - 1)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) >= abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je - 1) = -temp2r/temp
                       work(3*n + je - 1) = -temp2i/temp
                    else
                       work(2*n + je - 1) = one
                       work(3*n + je - 1) = zero
                       temp = acoef*s(je - 1,je)
                       work(2*n + je) = (bcoefr*p(je - 1,je - 1) - acoef*s(je - 1,je - 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je - 1,je - 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je - 1) &
                              ) + abs(work(3*n + je - 1)))
                    ! compute contribution from columns je and je-1
                    ! of a and b to the sums.
                    creala = acoef*work(2*n + je - 1)
                    cimaga = acoef*work(3*n + je - 1)
                    crealb = bcoefr*work(2*n + je - 1) - bcoefi*work(3*n + je - 1)
                    cimagb = bcoefi*work(2*n + je - 1) + bcoefr*work(3*n + je - 1)
                    cre2a = acoef*work(2*n + je)
                    cim2a = acoef*work(3*n + je)
                    cre2b = bcoefr*work(2*n + je) - bcoefi*work(3*n + je)
                    cim2b = bcoefi*work(2*n + je) + bcoefr*work(3*n + je)
                    do jr = 1,je - 2
                       work(2*n + jr) = -creala*s(jr,je - 1) + crealb*p(jr,je - 1) - cre2a*s(jr, &
                                 je) + cre2b*p(jr,je)
                       work(3*n + jr) = -cimaga*s(jr,je - 1) + cimagb*p(jr,je - 1) - cim2a*s(jr, &
                                 je) + cim2b*p(jr,je)
                    end do
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                 ! columnwise triangular solve of  (a a - b b)  x = 0
                 il2by2 = .false.
                 loop_370: do j = je - nw,1,-1
                    ! if a 2-by-2 block, is in position j-1:j, wait until
                    ! next iteration to process it (when it will be j:j+1)
                    if (.not. il2by2 .and. j > 1) then
                       if (s(j,j - 1) /= zero) then
                          il2by2 = .true.
                          cycle loop_370
                       end if
                    end if
                    bdiag(1) = p(j,j)
                    if (il2by2) then
                       na = 2
                       bdiag(2) = p(j + 1,j + 1)
                    else
                       na = 1
                    end if
                    ! compute x(j) (and x(j+1), if 2-by-2 block)
                    call la_xlaln2(.false.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),work(2*n + j),n,bcoefr,bcoefi,sum,2,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = 1,je
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                    end if
                    xmax = max(scale*xmax,temp)
                    do jw = 1,nw
                       do ja = 1,na
                          work((jw + 1)*n + j + ja - 1) = sum(ja,jw)
                       end do
                    end do
                    ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                    if (j > 1) then
                       ! check whether scaling is necessary for sum.
                       xscale = one/max(one,xmax)
                       temp = acoefa*work(j) + bcoefa*work(n + j)
                       if (il2by2) temp = max(temp,acoefa*work(j + 1) + bcoefa*work(n + j + 1))

                       temp = max(temp,acoefa,bcoefa)
                       if (temp > bignum*xscale) then
                          do jw = 0,nw - 1
                             do jr = 1,je
                                work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                             end do
                          end do
                          xmax = xmax*xscale
                       end if
                       ! compute the contributions of the off-diagonals of
                       ! column j (and j+1, if 2-by-2 block) of a and b to the
                       ! sums.
                       do ja = 1,na
                          if (ilcplx) then
                             creala = acoef*work(2*n + j + ja - 1)
                             cimaga = acoef*work(3*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1) - bcoefi*work(3*n + j + ja - 1)
                             cimagb = bcoefi*work(2*n + j + ja - 1) + bcoefr*work(3*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                                work(3*n + jr) = work(3*n + jr) - cimaga*s(jr,j + ja - 1) + cimagb*p( &
                                           jr,j + ja - 1)
                             end do
                          else
                             creala = acoef*work(2*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                             end do
                          end if
                       end do
                    end if
                    il2by2 = .false.
                 end do loop_370
                 ! copy eigenvector to vr, back transforming if
                 ! howmny='b'.
                 ieig = ieig - nw
                 if (ilback) then
                    do jw = 0,nw - 1
                       do jr = 1,n
                          work((jw + 4)*n + jr) = work((jw + 2)*n + 1)*vr(jr,1)
                       end do
                       ! a series of compiler directives to defeat
                       ! vectorization for the next loop
                       do jc = 2,je
                          do jr = 1,n
                             work((jw + 4)*n + jr) = work((jw + 4)*n + jr) + work((jw + 2)*n + jc) &
                                       *vr(jr,jc)
                          end do
                       end do
                    end do
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 4)*n + jr)
                       end do
                    end do
                    iend = n
                 else
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 2)*n + jr)
                       end do
                    end do
                    iend = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)) + abs(vr(j,ieig + 1)))
                    end do
                 else
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = 1,iend
                          vr(jr,ieig + jw) = xscale*vr(jr,ieig + jw)
                       end do
                    end do
                 end if
              end do loop_500
           end if
           return
     end subroutine la_xtgevc
#endif
#ifdef LA_WITH_QP
     !> QTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of real matrices (S,P), where S is a quasi-triangular matrix
     !> and P is upper triangular.  Matrix pairs of this type are produced by
     !> the generalized Schur factorization of a matrix pair (A,B):
     !> A = Q*S*Z**T,  B = Q*P*Z**T
     !> as computed by QGGHRD + QHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal blocks of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the orthogonal factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_qtgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(qp),intent(in) :: p(ldp,*),s(lds,*)
           real(qp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: safety = 1.0e+2_qp

           ! Local Scalars
           logical(lk) :: compl,compr,il2by2,ilabad,ilall,ilback,ilbbad,ilcomp,ilcplx, &
                     lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,iinfo,im,iside,j,ja,jc,je,jr,jw, &
                     na,nw
           real(qp) :: acoef,acoefa,anorm,ascale,bcoefa,bcoefi,bcoefr,big,bignum,bnorm, &
           bscale,cim2a,cim2b,cimaga,cimagb,cre2a,cre2b,creala,crealb,dmin,safmin, &
                     salfar,sbeta,scale,small,temp,temp2,temp2i,temp2r,ulp,xmax,xscale
           ! Local Arrays
           real(qp) :: bdiag(2),sum(2,2),sums(2,2),sump(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
              ilall = .true.
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors to be computed
           if (.not. ilall) then
              im = 0
              ilcplx = .false.
              loop_10: do j = 1,n
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_10
                 end if
                 if (j < n) then
                    if (s(j + 1,j) /= zero) ilcplx = .true.
                 end if
                 if (ilcplx) then
                    if (select(j) .or. select(j + 1)) im = im + 2
                 else
                    if (select(j)) im = im + 1
                 end if
              end do loop_10
           else
              im = n
           end if
           ! check 2-by-2 diagonal blocks of a, b
           ilabad = .false.
           ilbbad = .false.
           do j = 1,n - 1
              if (s(j + 1,j) /= zero) then
                 if (p(j,j) == zero .or. p(j + 1,j + 1) == zero .or. p(j,j + 1) /= zero) ilbbad = &
                           .true.
                 if (j < n - 1) then
                    if (s(j + 2,j + 1) /= zero) ilabad = .true.
                 end if
              end if
           end do
           if (ilabad) then
              info = -5
           else if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('QTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_qlamch('SAFE MINIMUM')
           big = one/safmin
           call la_qlabad(safmin,big)
           ulp = la_qlamch('EPSILON')*la_qlamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part (i.e., excluding all elements belonging to the diagonal
           ! blocks) of a and b to check for possible overflow in the
           ! triangular solver.
           anorm = abs(s(1,1))
           if (n > 1) anorm = anorm + abs(s(2,1))
           bnorm = abs(p(1,1))
           work(1) = zero
           work(n + 1) = zero
           do j = 2,n
              temp = zero
              temp2 = zero
              if (s(j,j - 1) == zero) then
                 iend = j - 1
              else
                 iend = j - 2
              end if
              do i = 1,iend
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              work(j) = temp
              work(n + j) = temp2
              do i = iend + 1,min(j + 1,n)
                 temp = temp + abs(s(i,j))
                 temp2 = temp2 + abs(p(i,j))
              end do
              anorm = max(anorm,temp)
              bnorm = max(bnorm,temp2)
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_220: do je = 1,n
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at.
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_220
                 end if
                 nw = 1
                 if (je < n) then
                    if (s(je + 1,je) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je + 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_220
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       ieig = ieig + 1
                       do jr = 1,n
                          vl(jr,ieig) = zero
                       end do
                       vl(ieig,ieig) = one
                       cycle loop_220
                    end if
                 end if
                 ! clear vector
                 do jr = 1,nw*n
                    work(2*n + jr) = zero
                 end do
                                                       ! t
                 ! compute coefficients in  ( a a - b b )  y = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                 else
                    ! complex eigenvalue
                    call la_qlag2(s(je,je),lds,p(je,je),ldp,safmin*safety,acoef, &
                              temp,bcoefr,temp2,bcoefi)
                    bcoefi = -bcoefi
                    if (bcoefi == zero) then
                       info = je
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    temp = acoef*s(je + 1,je)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) > abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je + 1) = -temp2r/temp
                       work(3*n + je + 1) = -temp2i/temp
                    else
                       work(2*n + je + 1) = one
                       work(3*n + je + 1) = zero
                       temp = acoef*s(je,je + 1)
                       work(2*n + je) = (bcoefr*p(je + 1,je + 1) - acoef*s(je + 1,je + 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je + 1,je + 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je + 1) &
                              ) + abs(work(3*n + je + 1)))
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                 ! t
                 ! triangular solve of  (a a - b b)  y = 0
                                         ! t
                 ! (rowwise in  (a a - b b) , or columnwise in (a a - b b) )
                 il2by2 = .false.
                 loop_160: do j = je + nw,n
                    if (il2by2) then
                       il2by2 = .false.
                       cycle loop_160
                    end if
                    na = 1
                    bdiag(1) = p(j,j)
                    if (j < n) then
                       if (s(j + 1,j) /= zero) then
                          il2by2 = .true.
                          bdiag(2) = p(j + 1,j + 1)
                          na = 2
                       end if
                    end if
                    ! check whether scaling is necessary for dot products
                    xscale = one/max(one,xmax)
                    temp = max(work(j),work(n + j),acoefa*work(j) + bcoefa*work(n + j))

                    if (il2by2) temp = max(temp,work(j + 1),work(n + j + 1),acoefa*work(j + 1) + &
                              bcoefa*work(n + j + 1))
                    if (temp > bignum*xscale) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = xmax*xscale
                    end if
                    ! compute dot products
                          ! j-1
                    ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                          ! k=je
                    ! to reduce the op count, this is done as
                    ! _        j-1                  _        j-1
                    ! a*conjg( sum  s(k,j)*x(k) ) - b*conjg( sum  p(k,j)*x(k) )
                             ! k=je                          k=je
                    ! which may cause underflow problems if a or b are close
                    ! to underflow.  (e.g., less than small.)
                    do jw = 1,nw
                       do ja = 1,na
                          sums(ja,jw) = zero
                          sump(ja,jw) = zero
                          do jr = je,j - 1
                             sums(ja,jw) = sums(ja,jw) + s(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                             sump(ja,jw) = sump(ja,jw) + p(jr,j + ja - 1)*work((jw + 1)*n + jr &
                                       )
                          end do
                       end do
                    end do
                    do ja = 1,na
                       if (ilcplx) then
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1) - bcoefi*sump( &
                                    ja,2)
                          sum(ja,2) = -acoef*sums(ja,2) + bcoefr*sump(ja,2) + bcoefi*sump( &
                                    ja,1)
                       else
                          sum(ja,1) = -acoef*sums(ja,1) + bcoefr*sump(ja,1)
                       end if
                    end do
                                        ! t
                    ! solve  ( a a - b b )  y = sum(,)
                    ! with scaling and perturbation of the denominator
                    call la_qlaln2(.true.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),sum,2,bcoefr,bcoefi,work(2*n + j),n,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = je,j - 1
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                       xmax = scale*xmax
                    end if
                    xmax = max(xmax,temp)
                 end do loop_160
                 ! copy eigenvector to vl, back transforming if
                 ! howmny='b'.
                 ieig = ieig + 1
                 if (ilback) then
                    do jw = 0,nw - 1
                       call la_qgemv('N',n,n + 1 - je,one,vl(1,je),ldvl,work((jw + 2)*n + &
                                 je),1,zero,work((jw + 4)*n + 1),1)
                    end do
                    call la_qlacpy(' ',n,nw,work(4*n + 1),n,vl(1,je),ldvl)
                    ibeg = 1
                 else
                    call la_qlacpy(' ',n,nw,work(2*n + 1),n,vl(1,ieig),ldvl)
                    ibeg = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)) + abs(vl(j,ieig + 1)))
                    end do
                 else
                    do j = ibeg,n
                       xmax = max(xmax,abs(vl(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = ibeg,n
                          vl(jr,ieig + jw) = xscale*vl(jr,ieig + jw)
                       end do
                    end do
                 end if
                 ieig = ieig + nw - 1
              end do loop_220
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              ilcplx = .false.
              loop_500: do je = n,1,-1
                 ! skip this iteration if (a) howmny='s' and select=.false., or
                 ! (b) this would be the second of a complex pair.
                 ! check for complex eigenvalue, so as to be sure of which
                 ! entry(-ies) of select to look at -- if complex, select(je)
                 ! or select(je-1).
                 ! if this is a complex pair, the 2-by-2 diagonal block
                 ! corresponding to the eigenvalue is in rows/columns je-1:je
                 if (ilcplx) then
                    ilcplx = .false.
                    cycle loop_500
                 end if
                 nw = 1
                 if (je > 1) then
                    if (s(je,je - 1) /= zero) then
                       ilcplx = .true.
                       nw = 2
                    end if
                 end if
                 if (ilall) then
                    ilcomp = .true.
                 else if (ilcplx) then
                    ilcomp = select(je) .or. select(je - 1)
                 else
                    ilcomp = select(je)
                 end if
                 if (.not. ilcomp) cycle loop_500
                 ! decide if (a) singular pencil, (b) real eigenvalue, or
                 ! (c) complex eigenvalue.
                 if (.not. ilcplx) then
                    if (abs(s(je,je)) <= safmin .and. abs(p(je,je)) <= safmin) then
                       ! singular matrix pencil -- unit eigenvector
                       ieig = ieig - 1
                       do jr = 1,n
                          vr(jr,ieig) = zero
                       end do
                       vr(ieig,ieig) = one
                       cycle loop_500
                    end if
                 end if
                 ! clear vector
                 do jw = 0,nw - 1
                    do jr = 1,n
                       work((jw + 2)*n + jr) = zero
                    end do
                 end do
                 ! compute coefficients in  ( a a - b b ) x = 0
                    ! a  is  acoef
                    ! b  is  bcoefr + i*bcoefi
                 if (.not. ilcplx) then
                    ! real eigenvalue
                    temp = one/max(abs(s(je,je))*ascale,abs(p(je,je))*bscale,safmin &
                              )
                    salfar = (temp*s(je,je))*ascale
                    sbeta = (temp*p(je,je))*bscale
                    acoef = sbeta*ascale
                    bcoefr = salfar*bscale
                    bcoefi = zero
                    ! scale to avoid underflow
                    scale = one
                    lsa = abs(sbeta) >= safmin .and. abs(acoef) < small
                    lsb = abs(salfar) >= safmin .and. abs(bcoefr) < small
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs(salfar))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoef),abs(bcoefr))) &
                                 )
                       if (lsa) then
                          acoef = ascale*(scale*sbeta)
                       else
                          acoef = scale*acoef
                       end if
                       if (lsb) then
                          bcoefr = bscale*(scale*salfar)
                       else
                          bcoefr = scale*bcoefr
                       end if
                    end if
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr)
                    ! first component is 1
                    work(2*n + je) = one
                    xmax = one
                    ! compute contribution from column je of a and b to sum
                    ! (see "further details", above.)
                    do jr = 1,je - 1
                       work(2*n + jr) = bcoefr*p(jr,je) - acoef*s(jr,je)
                    end do
                 else
                    ! complex eigenvalue
                    call la_qlag2(s(je - 1,je - 1),lds,p(je - 1,je - 1),ldp,safmin*safety, &
                              acoef,temp,bcoefr,temp2,bcoefi)
                    if (bcoefi == zero) then
                       info = je - 1
                       return
                    end if
                    ! scale to avoid over/underflow
                    acoefa = abs(acoef)
                    bcoefa = abs(bcoefr) + abs(bcoefi)
                    scale = one
                    if (acoefa*ulp < safmin .and. acoefa >= safmin) scale = (safmin/ulp)/ &
                              acoefa
                    if (bcoefa*ulp < safmin .and. bcoefa >= safmin) scale = max(scale, (safmin/ &
                              ulp)/bcoefa)
                    if (safmin*acoefa > ascale) scale = ascale/(safmin*acoefa)
                    if (safmin*bcoefa > bscale) scale = min(scale,bscale/(safmin*bcoefa))

                    if (scale /= one) then
                       acoef = scale*acoef
                       acoefa = abs(acoef)
                       bcoefr = scale*bcoefr
                       bcoefi = scale*bcoefi
                       bcoefa = abs(bcoefr) + abs(bcoefi)
                    end if
                    ! compute first two components of eigenvector
                    ! and contribution to sums
                    temp = acoef*s(je,je - 1)
                    temp2r = acoef*s(je,je) - bcoefr*p(je,je)
                    temp2i = -bcoefi*p(je,je)
                    if (abs(temp) >= abs(temp2r) + abs(temp2i)) then
                       work(2*n + je) = one
                       work(3*n + je) = zero
                       work(2*n + je - 1) = -temp2r/temp
                       work(3*n + je - 1) = -temp2i/temp
                    else
                       work(2*n + je - 1) = one
                       work(3*n + je - 1) = zero
                       temp = acoef*s(je - 1,je)
                       work(2*n + je) = (bcoefr*p(je - 1,je - 1) - acoef*s(je - 1,je - 1))/ &
                                 temp
                       work(3*n + je) = bcoefi*p(je - 1,je - 1)/temp
                    end if
                    xmax = max(abs(work(2*n + je)) + abs(work(3*n + je)),abs(work(2*n + je - 1) &
                              ) + abs(work(3*n + je - 1)))
                    ! compute contribution from columns je and je-1
                    ! of a and b to the sums.
                    creala = acoef*work(2*n + je - 1)
                    cimaga = acoef*work(3*n + je - 1)
                    crealb = bcoefr*work(2*n + je - 1) - bcoefi*work(3*n + je - 1)
                    cimagb = bcoefi*work(2*n + je - 1) + bcoefr*work(3*n + je - 1)
                    cre2a = acoef*work(2*n + je)
                    cim2a = acoef*work(3*n + je)
                    cre2b = bcoefr*work(2*n + je) - bcoefi*work(3*n + je)
                    cim2b = bcoefi*work(2*n + je) + bcoefr*work(3*n + je)
                    do jr = 1,je - 2
                       work(2*n + jr) = -creala*s(jr,je - 1) + crealb*p(jr,je - 1) - cre2a*s(jr, &
                                 je) + cre2b*p(jr,je)
                       work(3*n + jr) = -cimaga*s(jr,je - 1) + cimagb*p(jr,je - 1) - cim2a*s(jr, &
                                 je) + cim2b*p(jr,je)
                    end do
                 end if
                 dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                 ! columnwise triangular solve of  (a a - b b)  x = 0
                 il2by2 = .false.
                 loop_370: do j = je - nw,1,-1
                    ! if a 2-by-2 block, is in position j-1:j, wait until
                    ! next iteration to process it (when it will be j:j+1)
                    if (.not. il2by2 .and. j > 1) then
                       if (s(j,j - 1) /= zero) then
                          il2by2 = .true.
                          cycle loop_370
                       end if
                    end if
                    bdiag(1) = p(j,j)
                    if (il2by2) then
                       na = 2
                       bdiag(2) = p(j + 1,j + 1)
                    else
                       na = 1
                    end if
                    ! compute x(j) (and x(j+1), if 2-by-2 block)
                    call la_qlaln2(.false.,na,nw,dmin,acoef,s(j,j),lds,bdiag(1), &
                    bdiag(2),work(2*n + j),n,bcoefr,bcoefi,sum,2,scale,temp,iinfo)

                    if (scale < one) then
                       do jw = 0,nw - 1
                          do jr = 1,je
                             work((jw + 2)*n + jr) = scale*work((jw + 2)*n + jr)
                          end do
                       end do
                    end if
                    xmax = max(scale*xmax,temp)
                    do jw = 1,nw
                       do ja = 1,na
                          work((jw + 1)*n + j + ja - 1) = sum(ja,jw)
                       end do
                    end do
                    ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                    if (j > 1) then
                       ! check whether scaling is necessary for sum.
                       xscale = one/max(one,xmax)
                       temp = acoefa*work(j) + bcoefa*work(n + j)
                       if (il2by2) temp = max(temp,acoefa*work(j + 1) + bcoefa*work(n + j + 1))

                       temp = max(temp,acoefa,bcoefa)
                       if (temp > bignum*xscale) then
                          do jw = 0,nw - 1
                             do jr = 1,je
                                work((jw + 2)*n + jr) = xscale*work((jw + 2)*n + jr)
                             end do
                          end do
                          xmax = xmax*xscale
                       end if
                       ! compute the contributions of the off-diagonals of
                       ! column j (and j+1, if 2-by-2 block) of a and b to the
                       ! sums.
                       do ja = 1,na
                          if (ilcplx) then
                             creala = acoef*work(2*n + j + ja - 1)
                             cimaga = acoef*work(3*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1) - bcoefi*work(3*n + j + ja - 1)
                             cimagb = bcoefi*work(2*n + j + ja - 1) + bcoefr*work(3*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                                work(3*n + jr) = work(3*n + jr) - cimaga*s(jr,j + ja - 1) + cimagb*p( &
                                           jr,j + ja - 1)
                             end do
                          else
                             creala = acoef*work(2*n + j + ja - 1)
                             crealb = bcoefr*work(2*n + j + ja - 1)
                             do jr = 1,j - 1
                                work(2*n + jr) = work(2*n + jr) - creala*s(jr,j + ja - 1) + crealb*p( &
                                           jr,j + ja - 1)
                             end do
                          end if
                       end do
                    end if
                    il2by2 = .false.
                 end do loop_370
                 ! copy eigenvector to vr, back transforming if
                 ! howmny='b'.
                 ieig = ieig - nw
                 if (ilback) then
                    do jw = 0,nw - 1
                       do jr = 1,n
                          work((jw + 4)*n + jr) = work((jw + 2)*n + 1)*vr(jr,1)
                       end do
                       ! a series of compiler directives to defeat
                       ! vectorization for the next loop
                       do jc = 2,je
                          do jr = 1,n
                             work((jw + 4)*n + jr) = work((jw + 4)*n + jr) + work((jw + 2)*n + jc) &
                                       *vr(jr,jc)
                          end do
                       end do
                    end do
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 4)*n + jr)
                       end do
                    end do
                    iend = n
                 else
                    do jw = 0,nw - 1
                       do jr = 1,n
                          vr(jr,ieig + jw) = work((jw + 2)*n + jr)
                       end do
                    end do
                    iend = je
                 end if
                 ! scale eigenvector
                 xmax = zero
                 if (ilcplx) then
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)) + abs(vr(j,ieig + 1)))
                    end do
                 else
                    do j = 1,iend
                       xmax = max(xmax,abs(vr(j,ieig)))
                    end do
                 end if
                 if (xmax > safmin) then
                    xscale = one/xmax
                    do jw = 0,nw - 1
                       do jr = 1,iend
                          vr(jr,ieig + jw) = xscale*vr(jr,ieig + jw)
                       end do
                    end do
                 end if
              end do loop_500
           end if
           return
     end subroutine la_qtgevc
#endif

     !> STGEX2: swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
     !> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
     !> (A, B) by an orthogonal equivalence transformation.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,n1,n2, &
               work,lwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,lwork,n,n1,n2
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_scopy by calls to la_slaset, or by do
        ! loops. sven hammarling, 1/5/02.
           ! Parameters
           real(sp),parameter :: twenty = 2.0e+01_sp
           integer(ilp),parameter :: ldst = 4
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,idum,linfo,m
           real(sp) :: bqra21,brqa21,ddum,dnorma,dnormb,dscale,dsum,eps,f,g,sa,sb, &
                     scale,smlnum,thresha,threshb
           ! Local Arrays
           integer(ilp) :: iwork(ldst)
           real(sp) :: ai(2),ar(2),be(2),ir(ldst,ldst),ircop(ldst,ldst),li(ldst,ldst),licop( &
           ldst,ldst),s(ldst,ldst),scpy(ldst,ldst),t(ldst,ldst),taul(ldst),taur(ldst),tcpy( &
                     ldst,ldst)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1 .or. n1 <= 0 .or. n2 <= 0) return
           if (n1 > n .or. (j1 + n1) > n) return
           m = n1 + n2
           if (lwork < max(n*m,m*m*2)) then
              info = -16
              work(1) = max(n*m,m*m*2)
              return
           end if
           weak = .false.
           strong = .false.
           ! make a local copy of selected block
           call la_slaset('FULL',ldst,ldst,zero,zero,li,ldst)
           call la_slaset('FULL',ldst,ldst,zero,zero,ir,ldst)
           call la_slacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_slacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute threshold for testing acceptance of swapping.
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           dscale = zero
           dsum = one
           call la_slacpy('FULL',m,m,s,ldst,work,m)
           call la_slassq(m*m,work,1,dscale,dsum)
           dnorma = dscale*sqrt(dsum)
           dscale = zero
           dsum = one
           call la_slacpy('FULL',m,m,t,ldst,work,m)
           call la_slassq(m*m,work,1,dscale,dsum)
           dnormb = dscale*sqrt(dsum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*dnorma,smlnum)
           threshb = max(twenty*eps*dnormb,smlnum)
           if (m == 2) then
              ! case 1: swap 1-by-1 and 1-by-1 blocks.
              ! compute orthogonal ql and rq that swap 1-by-1 and 1-by-1 blocks
              ! using givens rotations and perform the swap tentatively.
              f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
              g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
              sa = abs(s(2,2))*abs(t(1,1))
              sb = abs(s(1,1))*abs(t(2,2))
              call la_slartg(f,g,ir(1,2),ir(1,1),ddum)
              ir(2,1) = -ir(1,2)
              ir(2,2) = ir(1,1)
              call la_srot(2,s(1,1),1,s(1,2),1,ir(1,1),ir(2,1))
              call la_srot(2,t(1,1),1,t(1,2),1,ir(1,1),ir(2,1))
              if (sa >= sb) then
                 call la_slartg(s(1,1),s(2,1),li(1,1),li(2,1),ddum)
              else
                 call la_slartg(t(1,1),t(2,1),li(1,1),li(2,1),ddum)
              end if
              call la_srot(2,s(1,1),ldst,s(2,1),ldst,li(1,1),li(2,1))

              call la_srot(2,t(1,1),ldst,t(2,1),ldst,li(1,1),li(2,1))

              li(2,2) = li(1,1)
              li(1,2) = -li(2,1)
              ! weak stability test: |s21| <= o(eps f-norm((a)))
                                 ! and  |t21| <= o(eps f-norm((b)))
              weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
              if (.not. weak) go to 70
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_slacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_sgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_sgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_slassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_slacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_sgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_sgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_slassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                     ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              call la_srot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_srot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_srot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,li(1,1),li(2,1 &
                        ))
              call la_srot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,li(1,1),li(2,1 &
                        ))
              ! set  n1-by-n2 (2,1) - blocks to zero.
              a(j1 + 1,j1) = zero
              b(j1 + 1,j1) = zero
              ! accumulate transformations into q and z if requested.
              if (wantz) call la_srot(n,z(1,j1),1,z(1,j1 + 1),1,ir(1,1),ir(2,1 &
                        ))
              if (wantq) call la_srot(n,q(1,j1),1,q(1,j1 + 1),1,li(1,1),li(2,1 &
                        ))
              ! exit with info = 0 if swap was successfully performed.
              return
           else
              ! case 2: swap 1-by-1 and 2-by-2 blocks, or 2-by-2
                      ! and 2-by-2 blocks.
              ! solve the generalized sylvester equation
                       ! s11 * r - l * s22 = scale * s12
                       ! t11 * r - l * t22 = scale * t12
              ! for r and l. solutions in li and ir.
              call la_slacpy('FULL',n1,n2,t(1,n1 + 1),ldst,li,ldst)
              call la_slacpy('FULL',n1,n2,s(1,n1 + 1),ldst,ir(n2 + 1,n1 + 1),ldst)

              call la_stgsy2('N',0,n1,n2,s,ldst,s(n1 + 1,n1 + 1),ldst,ir(n2 + 1,n1 + 1), &
               ldst,t,ldst,t(n1 + 1,n1 + 1),ldst,li,ldst,scale,dsum,dscale,iwork,idum, &
                         linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix ql:
                          ! ql**t * li = [ tl ]
                                       ! [ 0  ]
              ! where
                          ! li =  [      -l              ]
                                ! [ scale * identity(n2) ]
              do i = 1,n2
                 call la_sscal(n1,-one,li(1,i),1)
                 li(n1 + i,i) = scale
              end do
              call la_sgeqr2(m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_sorg2r(m,m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix rq:
                          ! ir * rq**t =   [ 0  tr],
               ! where ir = [ scale * identity(n1), r ]
              do i = 1,n1
                 ir(n2 + i,i) = scale
              end do
              call la_sgerq2(n1,m,ir(n2 + 1,1),ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_sorgr2(m,m,n1,ir,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              ! perform the swapping tentatively:
              call la_sgemm('T','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)
              call la_sgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,s,ldst)
              call la_sgemm('T','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)
              call la_sgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,t,ldst)
              call la_slacpy('F',m,m,s,ldst,scpy,ldst)
              call la_slacpy('F',m,m,t,ldst,tcpy,ldst)
              call la_slacpy('F',m,m,ir,ldst,ircop,ldst)
              call la_slacpy('F',m,m,li,ldst,licop,ldst)
              ! triangularize the b-part by an rq factorization.
              ! apply transformation (from left) to a-part, giving s.
              call la_sgerq2(m,m,t,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_sormr2('R','T',m,m,m,t,ldst,taur,s,ldst,work,linfo)
              if (linfo /= 0) go to 70
              call la_sormr2('L','N',m,m,m,t,ldst,taur,ir,ldst,work,linfo)
              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in brqa21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_slassq(n1,s(n2 + 1,i),1,dscale,dsum)
              end do
              brqa21 = dscale*sqrt(dsum)
              ! triangularize the b-part by a qr factorization.
              ! apply transformation (from right) to a-part, giving s.
              call la_sgeqr2(m,m,tcpy,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_sorm2r('L','T',m,m,m,tcpy,ldst,taul,scpy,ldst,work,info)

              call la_sorm2r('R','N',m,m,m,tcpy,ldst,taul,licop,ldst,work,info)

              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in bqra21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_slassq(n1,scpy(n2 + 1,i),1,dscale,dsum)
              end do
              bqra21 = dscale*sqrt(dsum)
              ! decide which method to use.
                ! weak stability test:
                   ! f-norm(s21) <= o(eps * f-norm((s)))
              if (bqra21 <= brqa21 .and. bqra21 <= thresha) then
                 call la_slacpy('F',m,m,scpy,ldst,s,ldst)
                 call la_slacpy('F',m,m,tcpy,ldst,t,ldst)
                 call la_slacpy('F',m,m,ircop,ldst,ir,ldst)
                 call la_slacpy('F',m,m,licop,ldst,li,ldst)
              else if (brqa21 >= thresha) then
                 go to 70
              end if
              ! set lower triangle of b-part to zero
              call la_slaset('LOWER',m - 1,m - 1,zero,zero,t(2,1),ldst)
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_slacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_sgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_sgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_slassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_slacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_sgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_sgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_slassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! if the swap is accepted ("weakly" and "strongly"), apply the
              ! transformations and set n1-by-n2 (2,1)-block to zero.
              call la_slaset('FULL',n1,n2,zero,zero,s(n2 + 1,1),ldst)
              ! copy back m-by-m diagonal block starting at index j1 of (a, b)
              call la_slacpy('F',m,m,s,ldst,a(j1,j1),lda)
              call la_slacpy('F',m,m,t,ldst,b(j1,j1),ldb)
              call la_slaset('FULL',ldst,ldst,zero,zero,t,ldst)
              ! standardize existing 2-by-2 blocks.
              call la_slaset('FULL',m,m,zero,zero,work,m)
              work(1) = one
              t(1,1) = one
              idum = lwork - m*m - 2
              if (n2 > 1) then
                 call la_slagv2(a(j1,j1),lda,b(j1,j1),ldb,ar,ai,be,work(1), &
                           work(2),t(1,1),t(2,1))
                 work(m + 1) = -work(2)
                 work(m + 2) = work(1)
                 t(n2,n2) = t(1,1)
                 t(1,2) = -t(2,1)
              end if
              work(m*m) = one
              t(m,m) = one
              if (n1 > 1) then
                 call la_slagv2(a(j1 + n2,j1 + n2),lda,b(j1 + n2,j1 + n2),ldb,taur,taul, &
                 work(m*m + 1),work(n2*m + n2 + 1),work(n2*m + n2 + 2),t(n2 + 1,n2 + 1),t(m,m - 1))

                 work(m*m) = work(n2*m + n2 + 1)
                 work(m*m - 1) = -work(n2*m + n2 + 2)
                 t(m,m) = t(n2 + 1,n2 + 1)
                 t(m - 1,m) = -t(m,m - 1)
              end if
              call la_sgemm('T','N',n2,n1,n2,one,work,m,a(j1,j1 + n2),lda,zero, &
                        work(m*m + 1),n2)
              call la_slacpy('FULL',n2,n1,work(m*m + 1),n2,a(j1,j1 + n2),lda)
              call la_sgemm('T','N',n2,n1,n2,one,work,m,b(j1,j1 + n2),ldb,zero, &
                        work(m*m + 1),n2)
              call la_slacpy('FULL',n2,n1,work(m*m + 1),n2,b(j1,j1 + n2),ldb)
              call la_sgemm('N','N',m,m,m,one,li,ldst,work,m,zero,work(m*m + 1),m &
                        )
              call la_slacpy('FULL',m,m,work(m*m + 1),m,li,ldst)
              call la_sgemm('N','N',n2,n1,n1,one,a(j1,j1 + n2),lda,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_slacpy('FULL',n2,n1,work,n2,a(j1,j1 + n2),lda)
              call la_sgemm('N','N',n2,n1,n1,one,b(j1,j1 + n2),ldb,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_slacpy('FULL',n2,n1,work,n2,b(j1,j1 + n2),ldb)
              call la_sgemm('T','N',m,m,m,one,ir,ldst,t,ldst,zero,work,m)
              call la_slacpy('FULL',m,m,work,m,ir,ldst)
              ! accumulate transformations into q and z if requested.
              if (wantq) then
                 call la_sgemm('N','N',n,m,m,one,q(1,j1),ldq,li,ldst,zero,work, &
                           n)
                 call la_slacpy('FULL',n,m,work,n,q(1,j1),ldq)
              end if
              if (wantz) then
                 call la_sgemm('N','N',n,m,m,one,z(1,j1),ldz,ir,ldst,zero,work, &
                           n)
                 call la_slacpy('FULL',n,m,work,n,z(1,j1),ldz)
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                      ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              i = j1 + m
              if (i <= n) then
                 call la_sgemm('T','N',m,n - i + 1,m,one,li,ldst,a(j1,i),lda,zero, &
                           work,m)
                 call la_slacpy('FULL',m,n - i + 1,work,m,a(j1,i),lda)
                 call la_sgemm('T','N',m,n - i + 1,m,one,li,ldst,b(j1,i),ldb,zero, &
                           work,m)
                 call la_slacpy('FULL',m,n - i + 1,work,m,b(j1,i),ldb)
              end if
              i = j1 - 1
              if (i > 0) then
                 call la_sgemm('N','N',i,m,m,one,a(1,j1),lda,ir,ldst,zero,work, &
                           i)
                 call la_slacpy('FULL',i,m,work,i,a(1,j1),lda)
                 call la_sgemm('N','N',i,m,m,one,b(1,j1),ldb,ir,ldst,zero,work, &
                           i)
                 call la_slacpy('FULL',i,m,work,i,b(1,j1),ldb)
              end if
              ! exit with info = 0 if swap was successfully performed.
              return
           end if
           ! exit with info = 1 if swap was rejected.
           70 continue
           info = 1
           return
     end subroutine la_stgex2
     !> DTGEX2: swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
     !> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
     !> (A, B) by an orthogonal equivalence transformation.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by DGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,n1,n2, &
               work,lwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,lwork,n,n1,n2
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_dcopy by calls to la_dlaset, or by do
        ! loops. sven hammarling, 1/5/02.
           ! Parameters
           real(dp),parameter :: twenty = 2.0e+01_dp
           integer(ilp),parameter :: ldst = 4
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,idum,linfo,m
           real(dp) :: bqra21,brqa21,ddum,dnorma,dnormb,dscale,dsum,eps,f,g,sa,sb, &
                     scale,smlnum,thresha,threshb
           ! Local Arrays
           integer(ilp) :: iwork(ldst)
           real(dp) :: ai(2),ar(2),be(2),ir(ldst,ldst),ircop(ldst,ldst),li(ldst,ldst),licop( &
           ldst,ldst),s(ldst,ldst),scpy(ldst,ldst),t(ldst,ldst),taul(ldst),taur(ldst),tcpy( &
                     ldst,ldst)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1 .or. n1 <= 0 .or. n2 <= 0) return
           if (n1 > n .or. (j1 + n1) > n) return
           m = n1 + n2
           if (lwork < max(1,n*m,m*m*2)) then
              info = -16
              work(1) = max(1,n*m,m*m*2)
              return
           end if
           weak = .false.
           strong = .false.
           ! make a local copy of selected block
           call la_dlaset('FULL',ldst,ldst,zero,zero,li,ldst)
           call la_dlaset('FULL',ldst,ldst,zero,zero,ir,ldst)
           call la_dlacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_dlacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute threshold for testing acceptance of swapping.
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           dscale = zero
           dsum = one
           call la_dlacpy('FULL',m,m,s,ldst,work,m)
           call la_dlassq(m*m,work,1,dscale,dsum)
           dnorma = dscale*sqrt(dsum)
           dscale = zero
           dsum = one
           call la_dlacpy('FULL',m,m,t,ldst,work,m)
           call la_dlassq(m*m,work,1,dscale,dsum)
           dnormb = dscale*sqrt(dsum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*dnorma,smlnum)
           threshb = max(twenty*eps*dnormb,smlnum)
           if (m == 2) then
              ! case 1: swap 1-by-1 and 1-by-1 blocks.
              ! compute orthogonal ql and rq that swap 1-by-1 and 1-by-1 blocks
              ! using givens rotations and perform the swap tentatively.
              f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
              g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
              sa = abs(s(2,2))*abs(t(1,1))
              sb = abs(s(1,1))*abs(t(2,2))
              call la_dlartg(f,g,ir(1,2),ir(1,1),ddum)
              ir(2,1) = -ir(1,2)
              ir(2,2) = ir(1,1)
              call la_drot(2,s(1,1),1,s(1,2),1,ir(1,1),ir(2,1))
              call la_drot(2,t(1,1),1,t(1,2),1,ir(1,1),ir(2,1))
              if (sa >= sb) then
                 call la_dlartg(s(1,1),s(2,1),li(1,1),li(2,1),ddum)
              else
                 call la_dlartg(t(1,1),t(2,1),li(1,1),li(2,1),ddum)
              end if
              call la_drot(2,s(1,1),ldst,s(2,1),ldst,li(1,1),li(2,1))

              call la_drot(2,t(1,1),ldst,t(2,1),ldst,li(1,1),li(2,1))

              li(2,2) = li(1,1)
              li(1,2) = -li(2,1)
              ! weak stability test: |s21| <= o(eps f-norm((a)))
                                 ! and  |t21| <= o(eps f-norm((b)))
              weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
              if (.not. weak) go to 70
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_dlacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_dgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_dgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_dlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_dlacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_dgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_dgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_dlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                     ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              call la_drot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_drot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_drot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,li(1,1),li(2,1 &
                        ))
              call la_drot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,li(1,1),li(2,1 &
                        ))
              ! set  n1-by-n2 (2,1) - blocks to zero.
              a(j1 + 1,j1) = zero
              b(j1 + 1,j1) = zero
              ! accumulate transformations into q and z if requested.
              if (wantz) call la_drot(n,z(1,j1),1,z(1,j1 + 1),1,ir(1,1),ir(2,1 &
                        ))
              if (wantq) call la_drot(n,q(1,j1),1,q(1,j1 + 1),1,li(1,1),li(2,1 &
                        ))
              ! exit with info = 0 if swap was successfully performed.
              return
           else
              ! case 2: swap 1-by-1 and 2-by-2 blocks, or 2-by-2
                      ! and 2-by-2 blocks.
              ! solve the generalized sylvester equation
                       ! s11 * r - l * s22 = scale * s12
                       ! t11 * r - l * t22 = scale * t12
              ! for r and l. solutions in li and ir.
              call la_dlacpy('FULL',n1,n2,t(1,n1 + 1),ldst,li,ldst)
              call la_dlacpy('FULL',n1,n2,s(1,n1 + 1),ldst,ir(n2 + 1,n1 + 1),ldst)

              call la_dtgsy2('N',0,n1,n2,s,ldst,s(n1 + 1,n1 + 1),ldst,ir(n2 + 1,n1 + 1), &
               ldst,t,ldst,t(n1 + 1,n1 + 1),ldst,li,ldst,scale,dsum,dscale,iwork,idum, &
                         linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix ql:
                          ! ql**t * li = [ tl ]
                                       ! [ 0  ]
              ! where
                          ! li =  [      -l              ]
                                ! [ scale * identity(n2) ]
              do i = 1,n2
                 call la_dscal(n1,-one,li(1,i),1)
                 li(n1 + i,i) = scale
              end do
              call la_dgeqr2(m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_dorg2r(m,m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix rq:
                          ! ir * rq**t =   [ 0  tr],
               ! where ir = [ scale * identity(n1), r ]
              do i = 1,n1
                 ir(n2 + i,i) = scale
              end do
              call la_dgerq2(n1,m,ir(n2 + 1,1),ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_dorgr2(m,m,n1,ir,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              ! perform the swapping tentatively:
              call la_dgemm('T','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)
              call la_dgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,s,ldst)
              call la_dgemm('T','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)
              call la_dgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,t,ldst)
              call la_dlacpy('F',m,m,s,ldst,scpy,ldst)
              call la_dlacpy('F',m,m,t,ldst,tcpy,ldst)
              call la_dlacpy('F',m,m,ir,ldst,ircop,ldst)
              call la_dlacpy('F',m,m,li,ldst,licop,ldst)
              ! triangularize the b-part by an rq factorization.
              ! apply transformation (from left) to a-part, giving s.
              call la_dgerq2(m,m,t,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_dormr2('R','T',m,m,m,t,ldst,taur,s,ldst,work,linfo)
              if (linfo /= 0) go to 70
              call la_dormr2('L','N',m,m,m,t,ldst,taur,ir,ldst,work,linfo)
              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in brqa21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_dlassq(n1,s(n2 + 1,i),1,dscale,dsum)
              end do
              brqa21 = dscale*sqrt(dsum)
              ! triangularize the b-part by a qr factorization.
              ! apply transformation (from right) to a-part, giving s.
              call la_dgeqr2(m,m,tcpy,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_dorm2r('L','T',m,m,m,tcpy,ldst,taul,scpy,ldst,work,info)

              call la_dorm2r('R','N',m,m,m,tcpy,ldst,taul,licop,ldst,work,info)

              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in bqra21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_dlassq(n1,scpy(n2 + 1,i),1,dscale,dsum)
              end do
              bqra21 = dscale*sqrt(dsum)
              ! decide which method to use.
                ! weak stability test:
                   ! f-norm(s21) <= o(eps * f-norm((s)))
              if (bqra21 <= brqa21 .and. bqra21 <= thresha) then
                 call la_dlacpy('F',m,m,scpy,ldst,s,ldst)
                 call la_dlacpy('F',m,m,tcpy,ldst,t,ldst)
                 call la_dlacpy('F',m,m,ircop,ldst,ir,ldst)
                 call la_dlacpy('F',m,m,licop,ldst,li,ldst)
              else if (brqa21 >= thresha) then
                 go to 70
              end if
              ! set lower triangle of b-part to zero
              call la_dlaset('LOWER',m - 1,m - 1,zero,zero,t(2,1),ldst)
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_dlacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_dgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_dgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_dlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_dlacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_dgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_dgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_dlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! if the swap is accepted ("weakly" and "strongly"), apply the
              ! transformations and set n1-by-n2 (2,1)-block to zero.
              call la_dlaset('FULL',n1,n2,zero,zero,s(n2 + 1,1),ldst)
              ! copy back m-by-m diagonal block starting at index j1 of (a, b)
              call la_dlacpy('F',m,m,s,ldst,a(j1,j1),lda)
              call la_dlacpy('F',m,m,t,ldst,b(j1,j1),ldb)
              call la_dlaset('FULL',ldst,ldst,zero,zero,t,ldst)
              ! standardize existing 2-by-2 blocks.
              call la_dlaset('FULL',m,m,zero,zero,work,m)
              work(1) = one
              t(1,1) = one
              idum = lwork - m*m - 2
              if (n2 > 1) then
                 call la_dlagv2(a(j1,j1),lda,b(j1,j1),ldb,ar,ai,be,work(1), &
                           work(2),t(1,1),t(2,1))
                 work(m + 1) = -work(2)
                 work(m + 2) = work(1)
                 t(n2,n2) = t(1,1)
                 t(1,2) = -t(2,1)
              end if
              work(m*m) = one
              t(m,m) = one
              if (n1 > 1) then
                 call la_dlagv2(a(j1 + n2,j1 + n2),lda,b(j1 + n2,j1 + n2),ldb,taur,taul, &
                 work(m*m + 1),work(n2*m + n2 + 1),work(n2*m + n2 + 2),t(n2 + 1,n2 + 1),t(m,m - 1))

                 work(m*m) = work(n2*m + n2 + 1)
                 work(m*m - 1) = -work(n2*m + n2 + 2)
                 t(m,m) = t(n2 + 1,n2 + 1)
                 t(m - 1,m) = -t(m,m - 1)
              end if
              call la_dgemm('T','N',n2,n1,n2,one,work,m,a(j1,j1 + n2),lda,zero, &
                        work(m*m + 1),n2)
              call la_dlacpy('FULL',n2,n1,work(m*m + 1),n2,a(j1,j1 + n2),lda)
              call la_dgemm('T','N',n2,n1,n2,one,work,m,b(j1,j1 + n2),ldb,zero, &
                        work(m*m + 1),n2)
              call la_dlacpy('FULL',n2,n1,work(m*m + 1),n2,b(j1,j1 + n2),ldb)
              call la_dgemm('N','N',m,m,m,one,li,ldst,work,m,zero,work(m*m + 1),m &
                        )
              call la_dlacpy('FULL',m,m,work(m*m + 1),m,li,ldst)
              call la_dgemm('N','N',n2,n1,n1,one,a(j1,j1 + n2),lda,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_dlacpy('FULL',n2,n1,work,n2,a(j1,j1 + n2),lda)
              call la_dgemm('N','N',n2,n1,n1,one,b(j1,j1 + n2),ldb,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_dlacpy('FULL',n2,n1,work,n2,b(j1,j1 + n2),ldb)
              call la_dgemm('T','N',m,m,m,one,ir,ldst,t,ldst,zero,work,m)
              call la_dlacpy('FULL',m,m,work,m,ir,ldst)
              ! accumulate transformations into q and z if requested.
              if (wantq) then
                 call la_dgemm('N','N',n,m,m,one,q(1,j1),ldq,li,ldst,zero,work, &
                           n)
                 call la_dlacpy('FULL',n,m,work,n,q(1,j1),ldq)
              end if
              if (wantz) then
                 call la_dgemm('N','N',n,m,m,one,z(1,j1),ldz,ir,ldst,zero,work, &
                           n)
                 call la_dlacpy('FULL',n,m,work,n,z(1,j1),ldz)
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                      ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              i = j1 + m
              if (i <= n) then
                 call la_dgemm('T','N',m,n - i + 1,m,one,li,ldst,a(j1,i),lda,zero, &
                           work,m)
                 call la_dlacpy('FULL',m,n - i + 1,work,m,a(j1,i),lda)
                 call la_dgemm('T','N',m,n - i + 1,m,one,li,ldst,b(j1,i),ldb,zero, &
                           work,m)
                 call la_dlacpy('FULL',m,n - i + 1,work,m,b(j1,i),ldb)
              end if
              i = j1 - 1
              if (i > 0) then
                 call la_dgemm('N','N',i,m,m,one,a(1,j1),lda,ir,ldst,zero,work, &
                           i)
                 call la_dlacpy('FULL',i,m,work,i,a(1,j1),lda)
                 call la_dgemm('N','N',i,m,m,one,b(1,j1),ldb,ir,ldst,zero,work, &
                           i)
                 call la_dlacpy('FULL',i,m,work,i,b(1,j1),ldb)
              end if
              ! exit with info = 0 if swap was successfully performed.
              return
           end if
           ! exit with info = 1 if swap was rejected.
           70 continue
           info = 1
           return
     end subroutine la_dtgex2
#ifdef LA_WITH_XDP
     !> XTGEX2: swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
     !> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
     !> (A, B) by an orthogonal equivalence transformation.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by XGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,n1,n2, &
               work,lwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,lwork,n,n1,n2
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_xcopy by calls to la_xlaset, or by do
        ! loops. sven hammarling, 1/5/02.
           ! Parameters
           real(xdp),parameter :: twenty = 2.0e+01_xdp
           integer(ilp),parameter :: ldst = 4
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,idum,linfo,m
           real(xdp) :: bqra21,brqa21,ddum,dnorma,dnormb,dscale,dsum,eps,f,g,sa,sb, &
                     scale,smlnum,thresha,threshb
           ! Local Arrays
           integer(ilp) :: iwork(ldst)
           real(xdp) :: ai(2),ar(2),be(2),ir(ldst,ldst),ircop(ldst,ldst),li(ldst,ldst),licop( &
           ldst,ldst),s(ldst,ldst),scpy(ldst,ldst),t(ldst,ldst),taul(ldst),taur(ldst),tcpy( &
                     ldst,ldst)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1 .or. n1 <= 0 .or. n2 <= 0) return
           if (n1 > n .or. (j1 + n1) > n) return
           m = n1 + n2
           if (lwork < max(1,n*m,m*m*2)) then
              info = -16
              work(1) = max(1,n*m,m*m*2)
              return
           end if
           weak = .false.
           strong = .false.
           ! make a local copy of selected block
           call la_xlaset('FULL',ldst,ldst,zero,zero,li,ldst)
           call la_xlaset('FULL',ldst,ldst,zero,zero,ir,ldst)
           call la_xlacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_xlacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute threshold for testing acceptance of swapping.
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           dscale = zero
           dsum = one
           call la_xlacpy('FULL',m,m,s,ldst,work,m)
           call la_xlassq(m*m,work,1,dscale,dsum)
           dnorma = dscale*sqrt(dsum)
           dscale = zero
           dsum = one
           call la_xlacpy('FULL',m,m,t,ldst,work,m)
           call la_xlassq(m*m,work,1,dscale,dsum)
           dnormb = dscale*sqrt(dsum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*dnorma,smlnum)
           threshb = max(twenty*eps*dnormb,smlnum)
           if (m == 2) then
              ! case 1: swap 1-by-1 and 1-by-1 blocks.
              ! compute orthogonal ql and rq that swap 1-by-1 and 1-by-1 blocks
              ! using givens rotations and perform the swap tentatively.
              f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
              g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
              sa = abs(s(2,2))*abs(t(1,1))
              sb = abs(s(1,1))*abs(t(2,2))
              call la_xlartg(f,g,ir(1,2),ir(1,1),ddum)
              ir(2,1) = -ir(1,2)
              ir(2,2) = ir(1,1)
              call la_xrot(2,s(1,1),1,s(1,2),1,ir(1,1),ir(2,1))
              call la_xrot(2,t(1,1),1,t(1,2),1,ir(1,1),ir(2,1))
              if (sa >= sb) then
                 call la_xlartg(s(1,1),s(2,1),li(1,1),li(2,1),ddum)
              else
                 call la_xlartg(t(1,1),t(2,1),li(1,1),li(2,1),ddum)
              end if
              call la_xrot(2,s(1,1),ldst,s(2,1),ldst,li(1,1),li(2,1))

              call la_xrot(2,t(1,1),ldst,t(2,1),ldst,li(1,1),li(2,1))

              li(2,2) = li(1,1)
              li(1,2) = -li(2,1)
              ! weak stability test: |s21| <= o(eps f-norm((a)))
                                 ! and  |t21| <= o(eps f-norm((b)))
              weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
              if (.not. weak) go to 70
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_xlacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_xgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_xgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_xlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_xlacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_xgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_xgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_xlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                     ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              call la_xrot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_xrot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_xrot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,li(1,1),li(2,1 &
                        ))
              call la_xrot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,li(1,1),li(2,1 &
                        ))
              ! set  n1-by-n2 (2,1) - blocks to zero.
              a(j1 + 1,j1) = zero
              b(j1 + 1,j1) = zero
              ! accumulate transformations into q and z if requested.
              if (wantz) call la_xrot(n,z(1,j1),1,z(1,j1 + 1),1,ir(1,1),ir(2,1 &
                        ))
              if (wantq) call la_xrot(n,q(1,j1),1,q(1,j1 + 1),1,li(1,1),li(2,1 &
                        ))
              ! exit with info = 0 if swap was successfully performed.
              return
           else
              ! case 2: swap 1-by-1 and 2-by-2 blocks, or 2-by-2
                      ! and 2-by-2 blocks.
              ! solve the generalized sylvester equation
                       ! s11 * r - l * s22 = scale * s12
                       ! t11 * r - l * t22 = scale * t12
              ! for r and l. solutions in li and ir.
              call la_xlacpy('FULL',n1,n2,t(1,n1 + 1),ldst,li,ldst)
              call la_xlacpy('FULL',n1,n2,s(1,n1 + 1),ldst,ir(n2 + 1,n1 + 1),ldst)

              call la_xtgsy2('N',0,n1,n2,s,ldst,s(n1 + 1,n1 + 1),ldst,ir(n2 + 1,n1 + 1), &
               ldst,t,ldst,t(n1 + 1,n1 + 1),ldst,li,ldst,scale,dsum,dscale,iwork,idum, &
                         linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix ql:
                          ! ql**t * li = [ tl ]
                                       ! [ 0  ]
              ! where
                          ! li =  [      -l              ]
                                ! [ scale * identity(n2) ]
              do i = 1,n2
                 call la_xscal(n1,-one,li(1,i),1)
                 li(n1 + i,i) = scale
              end do
              call la_xgeqr2(m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_xorg2r(m,m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix rq:
                          ! ir * rq**t =   [ 0  tr],
               ! where ir = [ scale * identity(n1), r ]
              do i = 1,n1
                 ir(n2 + i,i) = scale
              end do
              call la_xgerq2(n1,m,ir(n2 + 1,1),ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_xorgr2(m,m,n1,ir,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              ! perform the swapping tentatively:
              call la_xgemm('T','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)
              call la_xgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,s,ldst)
              call la_xgemm('T','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)
              call la_xgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,t,ldst)
              call la_xlacpy('F',m,m,s,ldst,scpy,ldst)
              call la_xlacpy('F',m,m,t,ldst,tcpy,ldst)
              call la_xlacpy('F',m,m,ir,ldst,ircop,ldst)
              call la_xlacpy('F',m,m,li,ldst,licop,ldst)
              ! triangularize the b-part by an rq factorization.
              ! apply transformation (from left) to a-part, giving s.
              call la_xgerq2(m,m,t,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_xormr2('R','T',m,m,m,t,ldst,taur,s,ldst,work,linfo)
              if (linfo /= 0) go to 70
              call la_xormr2('L','N',m,m,m,t,ldst,taur,ir,ldst,work,linfo)
              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in brqa21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_xlassq(n1,s(n2 + 1,i),1,dscale,dsum)
              end do
              brqa21 = dscale*sqrt(dsum)
              ! triangularize the b-part by a qr factorization.
              ! apply transformation (from right) to a-part, giving s.
              call la_xgeqr2(m,m,tcpy,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_xorm2r('L','T',m,m,m,tcpy,ldst,taul,scpy,ldst,work,info)

              call la_xorm2r('R','N',m,m,m,tcpy,ldst,taul,licop,ldst,work,info)

              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in bqra21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_xlassq(n1,scpy(n2 + 1,i),1,dscale,dsum)
              end do
              bqra21 = dscale*sqrt(dsum)
              ! decide which method to use.
                ! weak stability test:
                   ! f-norm(s21) <= o(eps * f-norm((s)))
              if (bqra21 <= brqa21 .and. bqra21 <= thresha) then
                 call la_xlacpy('F',m,m,scpy,ldst,s,ldst)
                 call la_xlacpy('F',m,m,tcpy,ldst,t,ldst)
                 call la_xlacpy('F',m,m,ircop,ldst,ir,ldst)
                 call la_xlacpy('F',m,m,licop,ldst,li,ldst)
              else if (brqa21 >= thresha) then
                 go to 70
              end if
              ! set lower triangle of b-part to zero
              call la_xlaset('LOWER',m - 1,m - 1,zero,zero,t(2,1),ldst)
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_xlacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_xgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_xgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_xlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_xlacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_xgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_xgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_xlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! if the swap is accepted ("weakly" and "strongly"), apply the
              ! transformations and set n1-by-n2 (2,1)-block to zero.
              call la_xlaset('FULL',n1,n2,zero,zero,s(n2 + 1,1),ldst)
              ! copy back m-by-m diagonal block starting at index j1 of (a, b)
              call la_xlacpy('F',m,m,s,ldst,a(j1,j1),lda)
              call la_xlacpy('F',m,m,t,ldst,b(j1,j1),ldb)
              call la_xlaset('FULL',ldst,ldst,zero,zero,t,ldst)
              ! standardize existing 2-by-2 blocks.
              call la_xlaset('FULL',m,m,zero,zero,work,m)
              work(1) = one
              t(1,1) = one
              idum = lwork - m*m - 2
              if (n2 > 1) then
                 call la_xlagv2(a(j1,j1),lda,b(j1,j1),ldb,ar,ai,be,work(1), &
                           work(2),t(1,1),t(2,1))
                 work(m + 1) = -work(2)
                 work(m + 2) = work(1)
                 t(n2,n2) = t(1,1)
                 t(1,2) = -t(2,1)
              end if
              work(m*m) = one
              t(m,m) = one
              if (n1 > 1) then
                 call la_xlagv2(a(j1 + n2,j1 + n2),lda,b(j1 + n2,j1 + n2),ldb,taur,taul, &
                 work(m*m + 1),work(n2*m + n2 + 1),work(n2*m + n2 + 2),t(n2 + 1,n2 + 1),t(m,m - 1))

                 work(m*m) = work(n2*m + n2 + 1)
                 work(m*m - 1) = -work(n2*m + n2 + 2)
                 t(m,m) = t(n2 + 1,n2 + 1)
                 t(m - 1,m) = -t(m,m - 1)
              end if
              call la_xgemm('T','N',n2,n1,n2,one,work,m,a(j1,j1 + n2),lda,zero, &
                        work(m*m + 1),n2)
              call la_xlacpy('FULL',n2,n1,work(m*m + 1),n2,a(j1,j1 + n2),lda)
              call la_xgemm('T','N',n2,n1,n2,one,work,m,b(j1,j1 + n2),ldb,zero, &
                        work(m*m + 1),n2)
              call la_xlacpy('FULL',n2,n1,work(m*m + 1),n2,b(j1,j1 + n2),ldb)
              call la_xgemm('N','N',m,m,m,one,li,ldst,work,m,zero,work(m*m + 1),m &
                        )
              call la_xlacpy('FULL',m,m,work(m*m + 1),m,li,ldst)
              call la_xgemm('N','N',n2,n1,n1,one,a(j1,j1 + n2),lda,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_xlacpy('FULL',n2,n1,work,n2,a(j1,j1 + n2),lda)
              call la_xgemm('N','N',n2,n1,n1,one,b(j1,j1 + n2),ldb,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_xlacpy('FULL',n2,n1,work,n2,b(j1,j1 + n2),ldb)
              call la_xgemm('T','N',m,m,m,one,ir,ldst,t,ldst,zero,work,m)
              call la_xlacpy('FULL',m,m,work,m,ir,ldst)
              ! accumulate transformations into q and z if requested.
              if (wantq) then
                 call la_xgemm('N','N',n,m,m,one,q(1,j1),ldq,li,ldst,zero,work, &
                           n)
                 call la_xlacpy('FULL',n,m,work,n,q(1,j1),ldq)
              end if
              if (wantz) then
                 call la_xgemm('N','N',n,m,m,one,z(1,j1),ldz,ir,ldst,zero,work, &
                           n)
                 call la_xlacpy('FULL',n,m,work,n,z(1,j1),ldz)
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                      ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              i = j1 + m
              if (i <= n) then
                 call la_xgemm('T','N',m,n - i + 1,m,one,li,ldst,a(j1,i),lda,zero, &
                           work,m)
                 call la_xlacpy('FULL',m,n - i + 1,work,m,a(j1,i),lda)
                 call la_xgemm('T','N',m,n - i + 1,m,one,li,ldst,b(j1,i),ldb,zero, &
                           work,m)
                 call la_xlacpy('FULL',m,n - i + 1,work,m,b(j1,i),ldb)
              end if
              i = j1 - 1
              if (i > 0) then
                 call la_xgemm('N','N',i,m,m,one,a(1,j1),lda,ir,ldst,zero,work, &
                           i)
                 call la_xlacpy('FULL',i,m,work,i,a(1,j1),lda)
                 call la_xgemm('N','N',i,m,m,one,b(1,j1),ldb,ir,ldst,zero,work, &
                           i)
                 call la_xlacpy('FULL',i,m,work,i,b(1,j1),ldb)
              end if
              ! exit with info = 0 if swap was successfully performed.
              return
           end if
           ! exit with info = 1 if swap was rejected.
           70 continue
           info = 1
           return
     end subroutine la_xtgex2
#endif
#ifdef LA_WITH_QP
     !> QTGEX2: swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
     !> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
     !> (A, B) by an orthogonal equivalence transformation.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by QGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,n1,n2, &
               work,lwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,lwork,n,n1,n2
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_qcopy by calls to la_qlaset, or by do
        ! loops. sven hammarling, 1/5/02.
           ! Parameters
           real(qp),parameter :: twenty = 2.0e+01_qp
           integer(ilp),parameter :: ldst = 4
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,idum,linfo,m
           real(qp) :: bqra21,brqa21,ddum,dnorma,dnormb,dscale,dsum,eps,f,g,sa,sb, &
                     scale,smlnum,thresha,threshb
           ! Local Arrays
           integer(ilp) :: iwork(ldst)
           real(qp) :: ai(2),ar(2),be(2),ir(ldst,ldst),ircop(ldst,ldst),li(ldst,ldst),licop( &
           ldst,ldst),s(ldst,ldst),scpy(ldst,ldst),t(ldst,ldst),taul(ldst),taur(ldst),tcpy( &
                     ldst,ldst)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1 .or. n1 <= 0 .or. n2 <= 0) return
           if (n1 > n .or. (j1 + n1) > n) return
           m = n1 + n2
           if (lwork < max(1,n*m,m*m*2)) then
              info = -16
              work(1) = max(1,n*m,m*m*2)
              return
           end if
           weak = .false.
           strong = .false.
           ! make a local copy of selected block
           call la_qlaset('FULL',ldst,ldst,zero,zero,li,ldst)
           call la_qlaset('FULL',ldst,ldst,zero,zero,ir,ldst)
           call la_qlacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_qlacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute threshold for testing acceptance of swapping.
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           dscale = zero
           dsum = one
           call la_qlacpy('FULL',m,m,s,ldst,work,m)
           call la_qlassq(m*m,work,1,dscale,dsum)
           dnorma = dscale*sqrt(dsum)
           dscale = zero
           dsum = one
           call la_qlacpy('FULL',m,m,t,ldst,work,m)
           call la_qlassq(m*m,work,1,dscale,dsum)
           dnormb = dscale*sqrt(dsum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*dnorma,smlnum)
           threshb = max(twenty*eps*dnormb,smlnum)
           if (m == 2) then
              ! case 1: swap 1-by-1 and 1-by-1 blocks.
              ! compute orthogonal ql and rq that swap 1-by-1 and 1-by-1 blocks
              ! using givens rotations and perform the swap tentatively.
              f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
              g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
              sa = abs(s(2,2))*abs(t(1,1))
              sb = abs(s(1,1))*abs(t(2,2))
              call la_qlartg(f,g,ir(1,2),ir(1,1),ddum)
              ir(2,1) = -ir(1,2)
              ir(2,2) = ir(1,1)
              call la_qrot(2,s(1,1),1,s(1,2),1,ir(1,1),ir(2,1))
              call la_qrot(2,t(1,1),1,t(1,2),1,ir(1,1),ir(2,1))
              if (sa >= sb) then
                 call la_qlartg(s(1,1),s(2,1),li(1,1),li(2,1),ddum)
              else
                 call la_qlartg(t(1,1),t(2,1),li(1,1),li(2,1),ddum)
              end if
              call la_qrot(2,s(1,1),ldst,s(2,1),ldst,li(1,1),li(2,1))

              call la_qrot(2,t(1,1),ldst,t(2,1),ldst,li(1,1),li(2,1))

              li(2,2) = li(1,1)
              li(1,2) = -li(2,1)
              ! weak stability test: |s21| <= o(eps f-norm((a)))
                                 ! and  |t21| <= o(eps f-norm((b)))
              weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
              if (.not. weak) go to 70
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_qlacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_qgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_qgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_qlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_qlacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_qgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_qgemm('N','T',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_qlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                     ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              call la_qrot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_qrot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,ir(1,1),ir(2,1))

              call la_qrot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,li(1,1),li(2,1 &
                        ))
              call la_qrot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,li(1,1),li(2,1 &
                        ))
              ! set  n1-by-n2 (2,1) - blocks to zero.
              a(j1 + 1,j1) = zero
              b(j1 + 1,j1) = zero
              ! accumulate transformations into q and z if requested.
              if (wantz) call la_qrot(n,z(1,j1),1,z(1,j1 + 1),1,ir(1,1),ir(2,1 &
                        ))
              if (wantq) call la_qrot(n,q(1,j1),1,q(1,j1 + 1),1,li(1,1),li(2,1 &
                        ))
              ! exit with info = 0 if swap was successfully performed.
              return
           else
              ! case 2: swap 1-by-1 and 2-by-2 blocks, or 2-by-2
                      ! and 2-by-2 blocks.
              ! solve the generalized sylvester equation
                       ! s11 * r - l * s22 = scale * s12
                       ! t11 * r - l * t22 = scale * t12
              ! for r and l. solutions in li and ir.
              call la_qlacpy('FULL',n1,n2,t(1,n1 + 1),ldst,li,ldst)
              call la_qlacpy('FULL',n1,n2,s(1,n1 + 1),ldst,ir(n2 + 1,n1 + 1),ldst)

              call la_qtgsy2('N',0,n1,n2,s,ldst,s(n1 + 1,n1 + 1),ldst,ir(n2 + 1,n1 + 1), &
               ldst,t,ldst,t(n1 + 1,n1 + 1),ldst,li,ldst,scale,dsum,dscale,iwork,idum, &
                         linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix ql:
                          ! ql**t * li = [ tl ]
                                       ! [ 0  ]
              ! where
                          ! li =  [      -l              ]
                                ! [ scale * identity(n2) ]
              do i = 1,n2
                 call la_qscal(n1,-one,li(1,i),1)
                 li(n1 + i,i) = scale
              end do
              call la_qgeqr2(m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_qorg2r(m,m,n2,li,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              ! compute orthogonal matrix rq:
                          ! ir * rq**t =   [ 0  tr],
               ! where ir = [ scale * identity(n1), r ]
              do i = 1,n1
                 ir(n2 + i,i) = scale
              end do
              call la_qgerq2(n1,m,ir(n2 + 1,1),ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_qorgr2(m,m,n1,ir,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              ! perform the swapping tentatively:
              call la_qgemm('T','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)
              call la_qgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,s,ldst)
              call la_qgemm('T','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)
              call la_qgemm('N','T',m,m,m,one,work,m,ir,ldst,zero,t,ldst)
              call la_qlacpy('F',m,m,s,ldst,scpy,ldst)
              call la_qlacpy('F',m,m,t,ldst,tcpy,ldst)
              call la_qlacpy('F',m,m,ir,ldst,ircop,ldst)
              call la_qlacpy('F',m,m,li,ldst,licop,ldst)
              ! triangularize the b-part by an rq factorization.
              ! apply transformation (from left) to a-part, giving s.
              call la_qgerq2(m,m,t,ldst,taur,work,linfo)
              if (linfo /= 0) go to 70
              call la_qormr2('R','T',m,m,m,t,ldst,taur,s,ldst,work,linfo)
              if (linfo /= 0) go to 70
              call la_qormr2('L','N',m,m,m,t,ldst,taur,ir,ldst,work,linfo)
              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in brqa21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_qlassq(n1,s(n2 + 1,i),1,dscale,dsum)
              end do
              brqa21 = dscale*sqrt(dsum)
              ! triangularize the b-part by a qr factorization.
              ! apply transformation (from right) to a-part, giving s.
              call la_qgeqr2(m,m,tcpy,ldst,taul,work,linfo)
              if (linfo /= 0) go to 70
              call la_qorm2r('L','T',m,m,m,tcpy,ldst,taul,scpy,ldst,work,info)

              call la_qorm2r('R','N',m,m,m,tcpy,ldst,taul,licop,ldst,work,info)

              if (linfo /= 0) go to 70
              ! compute f-norm(s21) in bqra21. (t21 is 0.)
              dscale = zero
              dsum = one
              do i = 1,n2
                 call la_qlassq(n1,scpy(n2 + 1,i),1,dscale,dsum)
              end do
              bqra21 = dscale*sqrt(dsum)
              ! decide which method to use.
                ! weak stability test:
                   ! f-norm(s21) <= o(eps * f-norm((s)))
              if (bqra21 <= brqa21 .and. bqra21 <= thresha) then
                 call la_qlacpy('F',m,m,scpy,ldst,s,ldst)
                 call la_qlacpy('F',m,m,tcpy,ldst,t,ldst)
                 call la_qlacpy('F',m,m,ircop,ldst,ir,ldst)
                 call la_qlacpy('F',m,m,licop,ldst,li,ldst)
              else if (brqa21 >= thresha) then
                 go to 70
              end if
              ! set lower triangle of b-part to zero
              call la_qlaset('LOWER',m - 1,m - 1,zero,zero,t(2,1),ldst)
              if (wands) then
                 ! strong stability test:
                     ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                     ! and
                     ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
                 call la_qlacpy('FULL',m,m,a(j1,j1),lda,work(m*m + 1),m)
                 call la_qgemm('N','N',m,m,m,one,li,ldst,s,ldst,zero,work,m)

                 call la_qgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_qlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sa = dscale*sqrt(dsum)
                 call la_qlacpy('FULL',m,m,b(j1,j1),ldb,work(m*m + 1),m)
                 call la_qgemm('N','N',m,m,m,one,li,ldst,t,ldst,zero,work,m)

                 call la_qgemm('N','N',m,m,m,-one,work,m,ir,ldst,one,work(m*m + 1), &
                            m)
                 dscale = zero
                 dsum = one
                 call la_qlassq(m*m,work(m*m + 1),1,dscale,dsum)
                 sb = dscale*sqrt(dsum)
                 strong = sa <= thresha .and. sb <= threshb
                 if (.not. strong) go to 70
              end if
              ! if the swap is accepted ("weakly" and "strongly"), apply the
              ! transformations and set n1-by-n2 (2,1)-block to zero.
              call la_qlaset('FULL',n1,n2,zero,zero,s(n2 + 1,1),ldst)
              ! copy back m-by-m diagonal block starting at index j1 of (a, b)
              call la_qlacpy('F',m,m,s,ldst,a(j1,j1),lda)
              call la_qlacpy('F',m,m,t,ldst,b(j1,j1),ldb)
              call la_qlaset('FULL',ldst,ldst,zero,zero,t,ldst)
              ! standardize existing 2-by-2 blocks.
              call la_qlaset('FULL',m,m,zero,zero,work,m)
              work(1) = one
              t(1,1) = one
              idum = lwork - m*m - 2
              if (n2 > 1) then
                 call la_qlagv2(a(j1,j1),lda,b(j1,j1),ldb,ar,ai,be,work(1), &
                           work(2),t(1,1),t(2,1))
                 work(m + 1) = -work(2)
                 work(m + 2) = work(1)
                 t(n2,n2) = t(1,1)
                 t(1,2) = -t(2,1)
              end if
              work(m*m) = one
              t(m,m) = one
              if (n1 > 1) then
                 call la_qlagv2(a(j1 + n2,j1 + n2),lda,b(j1 + n2,j1 + n2),ldb,taur,taul, &
                 work(m*m + 1),work(n2*m + n2 + 1),work(n2*m + n2 + 2),t(n2 + 1,n2 + 1),t(m,m - 1))

                 work(m*m) = work(n2*m + n2 + 1)
                 work(m*m - 1) = -work(n2*m + n2 + 2)
                 t(m,m) = t(n2 + 1,n2 + 1)
                 t(m - 1,m) = -t(m,m - 1)
              end if
              call la_qgemm('T','N',n2,n1,n2,one,work,m,a(j1,j1 + n2),lda,zero, &
                        work(m*m + 1),n2)
              call la_qlacpy('FULL',n2,n1,work(m*m + 1),n2,a(j1,j1 + n2),lda)
              call la_qgemm('T','N',n2,n1,n2,one,work,m,b(j1,j1 + n2),ldb,zero, &
                        work(m*m + 1),n2)
              call la_qlacpy('FULL',n2,n1,work(m*m + 1),n2,b(j1,j1 + n2),ldb)
              call la_qgemm('N','N',m,m,m,one,li,ldst,work,m,zero,work(m*m + 1),m &
                        )
              call la_qlacpy('FULL',m,m,work(m*m + 1),m,li,ldst)
              call la_qgemm('N','N',n2,n1,n1,one,a(j1,j1 + n2),lda,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_qlacpy('FULL',n2,n1,work,n2,a(j1,j1 + n2),lda)
              call la_qgemm('N','N',n2,n1,n1,one,b(j1,j1 + n2),ldb,t(n2 + 1,n2 + 1), &
                        ldst,zero,work,n2)
              call la_qlacpy('FULL',n2,n1,work,n2,b(j1,j1 + n2),ldb)
              call la_qgemm('T','N',m,m,m,one,ir,ldst,t,ldst,zero,work,m)
              call la_qlacpy('FULL',m,m,work,m,ir,ldst)
              ! accumulate transformations into q and z if requested.
              if (wantq) then
                 call la_qgemm('N','N',n,m,m,one,q(1,j1),ldq,li,ldst,zero,work, &
                           n)
                 call la_qlacpy('FULL',n,m,work,n,q(1,j1),ldq)
              end if
              if (wantz) then
                 call la_qgemm('N','N',n,m,m,one,z(1,j1),ldz,ir,ldst,zero,work, &
                           n)
                 call la_qlacpy('FULL',n,m,work,n,z(1,j1),ldz)
              end if
              ! update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
                      ! (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
              i = j1 + m
              if (i <= n) then
                 call la_qgemm('T','N',m,n - i + 1,m,one,li,ldst,a(j1,i),lda,zero, &
                           work,m)
                 call la_qlacpy('FULL',m,n - i + 1,work,m,a(j1,i),lda)
                 call la_qgemm('T','N',m,n - i + 1,m,one,li,ldst,b(j1,i),ldb,zero, &
                           work,m)
                 call la_qlacpy('FULL',m,n - i + 1,work,m,b(j1,i),ldb)
              end if
              i = j1 - 1
              if (i > 0) then
                 call la_qgemm('N','N',i,m,m,one,a(1,j1),lda,ir,ldst,zero,work, &
                           i)
                 call la_qlacpy('FULL',i,m,work,i,a(1,j1),lda)
                 call la_qgemm('N','N',i,m,m,one,b(1,j1),ldb,ir,ldst,zero,work, &
                           i)
                 call la_qlacpy('FULL',i,m,work,i,b(1,j1),ldb)
              end if
              ! exit with info = 0 if swap was successfully performed.
              return
           end if
           ! exit with info = 1 if swap was rejected.
           70 continue
           info = 1
           return
     end subroutine la_qtgex2
#endif

     !> STGEXC: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A,B) using an orthogonal equivalence transformation
     !> (A, B) = Q * (A, B) * Z**T,
     !> so that the diagonal block of (A, B) with row index IFST is moved
     !> to row ILST.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_stgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldq,ldz,lwork,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: here,lwmin,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info == 0) then
              if (n <= 1) then
                 lwmin = 1
              else
                 lwmin = 4*n + 16
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STGEXC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of the specified block and find out
           ! if it is 1-by-1 or 2-by-2.
           if (ifst > 1) then
              if (a(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (a(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out if it is 1-by-1 or 2-by-2.
           if (ilst > 1) then
              if (a(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (a(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst.
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (a(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,nbf, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (a(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here + 1,1, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,1, &
                              1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here + 1
                 else
                    ! recompute nbnext in case of 2-by-2 split.
                    if (a(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,nbnext,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2-by-2 block did split.
                       call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                       call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,nbf,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,1,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                              nbnext,1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here - 1
                 else
                   ! recompute nbnext in case of 2-by-2 split.
                    if (a(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - 1, &
                                  2,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2-by-2 block did split.
                       call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                       call la_stgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           work(1) = lwmin
           return
     end subroutine la_stgexc
     !> DTGEXC: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A,B) using an orthogonal equivalence transformation
     !> (A, B) = Q * (A, B) * Z**T,
     !> so that the diagonal block of (A, B) with row index IFST is moved
     !> to row ILST.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by DGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_dtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldq,ldz,lwork,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: here,lwmin,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info == 0) then
              if (n <= 1) then
                 lwmin = 1
              else
                 lwmin = 4*n + 16
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTGEXC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of the specified block and find out
           ! if it is 1-by-1 or 2-by-2.
           if (ifst > 1) then
              if (a(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (a(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out if it is 1-by-1 or 2-by-2.
           if (ilst > 1) then
              if (a(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (a(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst.
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (a(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,nbf, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (a(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here + 1,1, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,1, &
                              1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here + 1
                 else
                    ! recompute nbnext in case of 2-by-2 split.
                    if (a(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,nbnext,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2-by-2 block did split.
                       call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                       call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,nbf,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,1,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                              nbnext,1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here - 1
                 else
                   ! recompute nbnext in case of 2-by-2 split.
                    if (a(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - 1, &
                                  2,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2-by-2 block did split.
                       call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                       call la_dtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           work(1) = lwmin
           return
     end subroutine la_dtgexc
#ifdef LA_WITH_XDP
     !> XTGEXC: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A,B) using an orthogonal equivalence transformation
     !> (A, B) = Q * (A, B) * Z**T,
     !> so that the diagonal block of (A, B) with row index IFST is moved
     !> to row ILST.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by XGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_xtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldq,ldz,lwork,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: here,lwmin,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info == 0) then
              if (n <= 1) then
                 lwmin = 1
              else
                 lwmin = 4*n + 16
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTGEXC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of the specified block and find out
           ! if it is 1-by-1 or 2-by-2.
           if (ifst > 1) then
              if (a(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (a(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out if it is 1-by-1 or 2-by-2.
           if (ilst > 1) then
              if (a(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (a(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst.
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (a(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,nbf, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (a(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here + 1,1, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,1, &
                              1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here + 1
                 else
                    ! recompute nbnext in case of 2-by-2 split.
                    if (a(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,nbnext,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2-by-2 block did split.
                       call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                       call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,nbf,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,1,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                              nbnext,1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here - 1
                 else
                   ! recompute nbnext in case of 2-by-2 split.
                    if (a(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - 1, &
                                  2,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2-by-2 block did split.
                       call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                       call la_xtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           work(1) = lwmin
           return
     end subroutine la_xtgexc
#endif
#ifdef LA_WITH_QP
     !> QTGEXC: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A,B) using an orthogonal equivalence transformation
     !> (A, B) = Q * (A, B) * Z**T,
     !> so that the diagonal block of (A, B) with row index IFST is moved
     !> to row ILST.
     !> (A, B) must be in generalized real Schur canonical form (as returned
     !> by QGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
     !> diagonal blocks. B is upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
     !> Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T

     pure subroutine la_qtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldq,ldz,lwork,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: here,lwmin,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           lquery = (lwork == -1)
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info == 0) then
              if (n <= 1) then
                 lwmin = 1
              else
                 lwmin = 4*n + 16
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTGEXC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of the specified block and find out
           ! if it is 1-by-1 or 2-by-2.
           if (ifst > 1) then
              if (a(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (a(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out if it is 1-by-1 or 2-by-2.
           if (ilst > 1) then
              if (a(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (a(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst.
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (a(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,nbf, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (a(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here + 1,1, &
                           nbnext,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,1, &
                              1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here + 1
                 else
                    ! recompute nbnext in case of 2-by-2 split.
                    if (a(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,nbnext,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2-by-2 block did split.
                       call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                       call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 1
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap with next one below.
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1-by-1 or 2-by-2.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,nbf,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2-by-2 block breaks into two 1-by-1 blocks.
                 if (nbf == 2) then
                    if (a(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1-by-1 blocks, each of which
                 ! must be swapped individually.
                 nbnext = 1
                 if (here >= 3) then
                    if (a(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - nbnext, &
                           nbnext,1,work,lwork,info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1-by-1 blocks.
                    call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                              nbnext,1,work,lwork,info)
                    if (info /= 0) then
                       ilst = here
                       return
                    end if
                    here = here - 1
                 else
                   ! recompute nbnext in case of 2-by-2 split.
                    if (a(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2-by-2 block did not split.
                       call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here - 1, &
                                  2,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2-by-2 block did split.
                       call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                       call la_qtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here, &
                                 1,1,work,lwork,info)
                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 1
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           work(1) = lwmin
           return
     end subroutine la_qtgexc
#endif

     !> STGSEN: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A, B) (in terms of an orthonormal equivalence trans-
     !> formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the upper quasi-triangular
     !> matrix A and the upper triangular B. The leading columns of Q and
     !> Z form orthonormal bases of the corresponding left and right eigen-
     !> spaces (deflating subspaces). (A, B) must be in generalized real
     !> Schur canonical form (as returned by SGGES), i.e. A is block upper
     !> triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
     !> triangular.
     !> STGSEN also computes the generalized eigenvalues
     !> w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, STGSEN computes the estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_stgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alphar,alphai, &
               beta,q,ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),dif(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,kk,ks,liwmin,lwmin,mn2,n1,n2
           real(sp) :: dscale,dsum,eps,rdscal,smlnum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -14
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('STGSEN',-info)
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           pair = .false.
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) == zero) then
                       if (select(k)) m = m + 1
                    else
                       pair = .true.
                       if (select(k) .or. select(k + 1)) m = m + 2
                    end if
                 else
                    if (select(n)) m = m + 1
                 end if
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,4*n + 16,2*m*(n - m))
              liwmin = max(1,n + 6)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*n + 16,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 6)
           else
              lwmin = max(1,4*n + 16)
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -22
           else if (liwork < liwmin .and. .not. lquery) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('STGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_slassq(n,a(1,i),1,dscale,dsum)
                    call la_slassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 60
           end if
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           pair = .false.
           loop_30: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ! perform the reordering of diagonal blocks in (a, b)
                    ! by orthogonal transformation matrices and update
                    ! q and z accordingly (if requested):
                    kk = k
                    if (k /= ks) call la_stgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz, &
                               kk,ks,work,lwork,ierr)
                    if (ierr > 0) then
                       ! swap is rejected: exit.
                       info = 1
                       if (wantp) then
                          pl = zero
                          pr = zero
                       end if
                       if (wantd) then
                          dif(1) = zero
                          dif(2) = zero
                       end if
                       go to 60
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_30
           if (wantp) then
              ! solve generalized sylvester equation for r and l
              ! and compute pl and pr.
              n1 = m
              n2 = n - m
              i = n1 + 1
              ijb = 0
              call la_slacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_slacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              call la_stgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto left
              ! and right eigenspaces.
              rdscal = zero
              dsum = one
              call la_slassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_slassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates of difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu-estimate.
                 call la_stgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl-estimate.
                 call la_stgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_slacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_slacn2(mn2,work(mn2 + 1),work,iwork,dif(1),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_stgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_stgsyl('T',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_slacn2(mn2,work(mn2 + 1),work,iwork,dif(2),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_stgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_stgsyl('T',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           60 continue
           ! compute generalized eigenvalues of reordered pair (a, b) and
           ! normalize the generalized schur form.
           pair = .false.
           loop_70: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                    end if
                 end if
                 if (pair) then
                   ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_slag2(work,2,work(5),2,smlnum*eps,beta(k),beta(k + 1), &
                              alphar(k),alphar(k + 1),alphai(k))
                    alphai(k + 1) = -alphai(k)
                 else
                    if (sign(one,b(k,k)) < zero) then
                       ! if b(k,k) is negative, make it positive
                       do i = 1,n
                          a(k,i) = -a(k,i)
                          b(k,i) = -b(k,i)
                          if (wantq) q(i,k) = -q(i,k)
                       end do
                    end if
                    alphar(k) = a(k,k)
                    alphai(k) = zero
                    beta(k) = b(k,k)
                 end if
              end if
           end do loop_70
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_stgsen
     !> DTGSEN: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A, B) (in terms of an orthonormal equivalence trans-
     !> formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the upper quasi-triangular
     !> matrix A and the upper triangular B. The leading columns of Q and
     !> Z form orthonormal bases of the corresponding left and right eigen-
     !> spaces (deflating subspaces). (A, B) must be in generalized real
     !> Schur canonical form (as returned by DGGES), i.e. A is block upper
     !> triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
     !> triangular.
     !> DTGSEN also computes the generalized eigenvalues
     !> w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, DTGSEN computes the estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_dtgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alphar,alphai, &
               beta,q,ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),dif(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,kk,ks,liwmin,lwmin,mn2,n1,n2
           real(dp) :: dscale,dsum,eps,rdscal,smlnum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -14
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('DTGSEN',-info)
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           pair = .false.
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) == zero) then
                       if (select(k)) m = m + 1
                    else
                       pair = .true.
                       if (select(k) .or. select(k + 1)) m = m + 2
                    end if
                 else
                    if (select(n)) m = m + 1
                 end if
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,4*n + 16,2*m*(n - m))
              liwmin = max(1,n + 6)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*n + 16,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 6)
           else
              lwmin = max(1,4*n + 16)
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -22
           else if (liwork < liwmin .and. .not. lquery) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('DTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_dlassq(n,a(1,i),1,dscale,dsum)
                    call la_dlassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 60
           end if
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           pair = .false.
           loop_30: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ! perform the reordering of diagonal blocks in (a, b)
                    ! by orthogonal transformation matrices and update
                    ! q and z accordingly (if requested):
                    kk = k
                    if (k /= ks) call la_dtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz, &
                               kk,ks,work,lwork,ierr)
                    if (ierr > 0) then
                       ! swap is rejected: exit.
                       info = 1
                       if (wantp) then
                          pl = zero
                          pr = zero
                       end if
                       if (wantd) then
                          dif(1) = zero
                          dif(2) = zero
                       end if
                       go to 60
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_30
           if (wantp) then
              ! solve generalized sylvester equation for r and l
              ! and compute pl and pr.
              n1 = m
              n2 = n - m
              i = n1 + 1
              ijb = 0
              call la_dlacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_dlacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              call la_dtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto left
              ! and right eigenspaces.
              rdscal = zero
              dsum = one
              call la_dlassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_dlassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates of difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu-estimate.
                 call la_dtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl-estimate.
                 call la_dtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_dlacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_dlacn2(mn2,work(mn2 + 1),work,iwork,dif(1),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_dtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_dtgsyl('T',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_dlacn2(mn2,work(mn2 + 1),work,iwork,dif(2),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_dtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_dtgsyl('T',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           60 continue
           ! compute generalized eigenvalues of reordered pair (a, b) and
           ! normalize the generalized schur form.
           pair = .false.
           loop_80: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                    end if
                 end if
                 if (pair) then
                   ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_dlag2(work,2,work(5),2,smlnum*eps,beta(k),beta(k + 1), &
                              alphar(k),alphar(k + 1),alphai(k))
                    alphai(k + 1) = -alphai(k)
                 else
                    if (sign(one,b(k,k)) < zero) then
                       ! if b(k,k) is negative, make it positive
                       do i = 1,n
                          a(k,i) = -a(k,i)
                          b(k,i) = -b(k,i)
                          if (wantq) q(i,k) = -q(i,k)
                       end do
                    end if
                    alphar(k) = a(k,k)
                    alphai(k) = zero
                    beta(k) = b(k,k)
                 end if
              end if
           end do loop_80
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_dtgsen
#ifdef LA_WITH_XDP
     !> XTGSEN: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A, B) (in terms of an orthonormal equivalence trans-
     !> formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the upper quasi-triangular
     !> matrix A and the upper triangular B. The leading columns of Q and
     !> Z form orthonormal bases of the corresponding left and right eigen-
     !> spaces (deflating subspaces). (A, B) must be in generalized real
     !> Schur canonical form (as returned by XGGES), i.e. A is block upper
     !> triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
     !> triangular.
     !> XTGSEN also computes the generalized eigenvalues
     !> w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, XTGSEN computes the estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_xtgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alphar,alphai, &
               beta,q,ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(xdp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(xdp),intent(out) :: alphai(*),alphar(*),beta(*),dif(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,kk,ks,liwmin,lwmin,mn2,n1,n2
           real(xdp) :: dscale,dsum,eps,rdscal,smlnum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -14
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('XTGSEN',-info)
              return
           end if
           ! get machine constants
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           pair = .false.
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) == zero) then
                       if (select(k)) m = m + 1
                    else
                       pair = .true.
                       if (select(k) .or. select(k + 1)) m = m + 2
                    end if
                 else
                    if (select(n)) m = m + 1
                 end if
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,4*n + 16,2*m*(n - m))
              liwmin = max(1,n + 6)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*n + 16,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 6)
           else
              lwmin = max(1,4*n + 16)
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -22
           else if (liwork < liwmin .and. .not. lquery) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('XTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_xlassq(n,a(1,i),1,dscale,dsum)
                    call la_xlassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 60
           end if
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           pair = .false.
           loop_30: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ! perform the reordering of diagonal blocks in (a, b)
                    ! by orthogonal transformation matrices and update
                    ! q and z accordingly (if requested):
                    kk = k
                    if (k /= ks) call la_xtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz, &
                               kk,ks,work,lwork,ierr)
                    if (ierr > 0) then
                       ! swap is rejected: exit.
                       info = 1
                       if (wantp) then
                          pl = zero
                          pr = zero
                       end if
                       if (wantd) then
                          dif(1) = zero
                          dif(2) = zero
                       end if
                       go to 60
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_30
           if (wantp) then
              ! solve generalized sylvester equation for r and l
              ! and compute pl and pr.
              n1 = m
              n2 = n - m
              i = n1 + 1
              ijb = 0
              call la_xlacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_xlacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              call la_xtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto left
              ! and right eigenspaces.
              rdscal = zero
              dsum = one
              call la_xlassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_xlassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates of difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu-estimate.
                 call la_xtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl-estimate.
                 call la_xtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_xlacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_xlacn2(mn2,work(mn2 + 1),work,iwork,dif(1),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_xtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_xtgsyl('T',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_xlacn2(mn2,work(mn2 + 1),work,iwork,dif(2),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_xtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_xtgsyl('T',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           60 continue
           ! compute generalized eigenvalues of reordered pair (a, b) and
           ! normalize the generalized schur form.
           pair = .false.
           loop_80: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                    end if
                 end if
                 if (pair) then
                   ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_xlag2(work,2,work(5),2,smlnum*eps,beta(k),beta(k + 1), &
                              alphar(k),alphar(k + 1),alphai(k))
                    alphai(k + 1) = -alphai(k)
                 else
                    if (sign(one,b(k,k)) < zero) then
                       ! if b(k,k) is negative, make it positive
                       do i = 1,n
                          a(k,i) = -a(k,i)
                          b(k,i) = -b(k,i)
                          if (wantq) q(i,k) = -q(i,k)
                       end do
                    end if
                    alphar(k) = a(k,k)
                    alphai(k) = zero
                    beta(k) = b(k,k)
                 end if
              end if
           end do loop_80
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_xtgsen
#endif
#ifdef LA_WITH_QP
     !> QTGSEN: reorders the generalized real Schur decomposition of a real
     !> matrix pair (A, B) (in terms of an orthonormal equivalence trans-
     !> formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the upper quasi-triangular
     !> matrix A and the upper triangular B. The leading columns of Q and
     !> Z form orthonormal bases of the corresponding left and right eigen-
     !> spaces (deflating subspaces). (A, B) must be in generalized real
     !> Schur canonical form (as returned by QGGES), i.e. A is block upper
     !> triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
     !> triangular.
     !> QTGSEN also computes the generalized eigenvalues
     !> w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, QTGSEN computes the estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_qtgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alphar,alphai, &
               beta,q,ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),dif(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,kk,ks,liwmin,lwmin,mn2,n1,n2
           real(qp) :: dscale,dsum,eps,rdscal,smlnum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -14
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('QTGSEN',-info)
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           pair = .false.
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) == zero) then
                       if (select(k)) m = m + 1
                    else
                       pair = .true.
                       if (select(k) .or. select(k + 1)) m = m + 2
                    end if
                 else
                    if (select(n)) m = m + 1
                 end if
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,4*n + 16,2*m*(n - m))
              liwmin = max(1,n + 6)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*n + 16,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 6)
           else
              lwmin = max(1,4*n + 16)
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -22
           else if (liwork < liwmin .and. .not. lquery) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('QTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_qlassq(n,a(1,i),1,dscale,dsum)
                    call la_qlassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 60
           end if
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           pair = .false.
           loop_30: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ! perform the reordering of diagonal blocks in (a, b)
                    ! by orthogonal transformation matrices and update
                    ! q and z accordingly (if requested):
                    kk = k
                    if (k /= ks) call la_qtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz, &
                               kk,ks,work,lwork,ierr)
                    if (ierr > 0) then
                       ! swap is rejected: exit.
                       info = 1
                       if (wantp) then
                          pl = zero
                          pr = zero
                       end if
                       if (wantd) then
                          dif(1) = zero
                          dif(2) = zero
                       end if
                       go to 60
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_30
           if (wantp) then
              ! solve generalized sylvester equation for r and l
              ! and compute pl and pr.
              n1 = m
              n2 = n - m
              i = n1 + 1
              ijb = 0
              call la_qlacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_qlacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              call la_qtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto left
              ! and right eigenspaces.
              rdscal = zero
              dsum = one
              call la_qlassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_qlassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates of difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu-estimate.
                 call la_qtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl-estimate.
                 call la_qtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_qlacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_qlacn2(mn2,work(mn2 + 1),work,iwork,dif(1),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_qtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_qtgsyl('T',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_qlacn2(mn2,work(mn2 + 1),work,iwork,dif(2),kase,isave)

                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation.
                       call la_qtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_qtgsyl('T',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(2*n1*n2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           60 continue
           ! compute generalized eigenvalues of reordered pair (a, b) and
           ! normalize the generalized schur form.
           pair = .false.
           loop_80: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 if (k < n) then
                    if (a(k + 1,k) /= zero) then
                       pair = .true.
                    end if
                 end if
                 if (pair) then
                   ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_qlag2(work,2,work(5),2,smlnum*eps,beta(k),beta(k + 1), &
                              alphar(k),alphar(k + 1),alphai(k))
                    alphai(k + 1) = -alphai(k)
                 else
                    if (sign(one,b(k,k)) < zero) then
                       ! if b(k,k) is negative, make it positive
                       do i = 1,n
                          a(k,i) = -a(k,i)
                          b(k,i) = -b(k,i)
                          if (wantq) q(i,k) = -q(i,k)
                       end do
                    end if
                    alphar(k) = a(k,k)
                    alphai(k) = zero
                    beta(k) = b(k,k)
                 end if
              end if
           end do loop_80
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_qtgsen
#endif

     !> STGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B) in
     !> generalized real Schur canonical form (or of any matrix pair
     !> (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where
     !> Z**T denotes the transpose of Z.
     !> (A, B) must be in generalized real Schur form (as returned by SGGES),
     !> i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
     !> blocks. B is upper triangular.

     pure subroutine la_stgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_sp,only:zero,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           real(sp),intent(out) :: dif(*),s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: difdri = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,iz,k,ks,lwmin,n1,n2
           real(sp) :: alphai,alphar,alprqt,beta,c1,c2,cond,eps,lnrm,rnrm,root1,root2, &
                     scale,smlnum,tmpii,tmpir,tmpri,tmprr,uhav,uhavi,uhbv,uhbvi
           ! Local Arrays
           real(sp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 pair = .false.
                 do k = 1,n
                    if (pair) then
                       pair = .false.
                    else
                       if (k < n) then
                          if (a(k + 1,k) == zero) then
                             if (select(k)) m = m + 1
                          else
                             pair = .true.
                             if (select(k) .or. select(k + 1)) m = m + 2
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*(n + 2) + 16
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              ! determine whether a(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_20
              else
                 if (k < n) pair = a(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_20
                 else
                    if (.not. select(k)) cycle loop_20
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (pair) then
                    ! complex eigenvalue pair.
                    rnrm = la_slapy2(la_snrm2(n,vr(1,ks),1),la_snrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_slapy2(la_snrm2(n,vl(1,ks),1),la_snrm2(n,vl( &
                              1,ks + 1),1))
                    call la_sgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    tmprr = la_sdot(n,work,1,vl(1,ks),1)
                    tmpri = la_sdot(n,work,1,vl(1,ks + 1),1)
                    call la_sgemv('N',n,n,one,a,lda,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_sdot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_sdot(n,work,1,vl(1,ks),1)
                    uhav = tmprr + tmpii
                    uhavi = tmpir - tmpri
                    call la_sgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    tmprr = la_sdot(n,work,1,vl(1,ks),1)
                    tmpri = la_sdot(n,work,1,vl(1,ks + 1),1)
                    call la_sgemv('N',n,n,one,b,ldb,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_sdot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_sdot(n,work,1,vl(1,ks),1)
                    uhbv = tmprr + tmpii
                    uhbvi = tmpir - tmpri
                    uhav = la_slapy2(uhav,uhavi)
                    uhbv = la_slapy2(uhbv,uhbvi)
                    cond = la_slapy2(uhav,uhbv)
                    s(ks) = cond/(rnrm*lnrm)
                    s(ks + 1) = s(ks)
                 else
                    ! real eigenvalue.
                    rnrm = la_snrm2(n,vr(1,ks),1)
                    lnrm = la_snrm2(n,vl(1,ks),1)
                    call la_sgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    uhav = la_sdot(n,work,1,vl(1,ks),1)
                    call la_sgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    uhbv = la_sdot(n,work,1,vl(1,ks),1)
                    cond = la_slapy2(uhav,uhbv)
                    if (cond == zero) then
                       s(ks) = -one
                    else
                       s(ks) = cond/(rnrm*lnrm)
                    end if
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_slapy2(a(1,1),b(1,1))
                    cycle loop_20
                 end if
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvectors.
                 if (pair) then
                    ! copy the  2-by 2 pencil beginning at (a(k,k), b(k, k)).
                    ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_slag2(work,2,work(5),2,smlnum*eps,beta,dummy1(1), &
                              alphar,dummy(1),alphai)
                    alprqt = one
                    c1 = two*(alphar*alphar + alphai*alphai + beta*beta)
                    c2 = four*beta*beta*alphai*alphai
                    root1 = c1 + sqrt(c1*c1 - 4.0_sp*c2)
                    root2 = c2/root1
                    root1 = root1/two
                    cond = min(sqrt(root1),sqrt(root2))
                 end if
                 ! copy the matrix (a, b) to the array work and swap the
                 ! diagonal block beginning at a(k,k) to the (1,1) position.
                 call la_slacpy('FULL',n,n,a,lda,work,n)
                 call la_slacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                 ifst = k
                 ilst = 1
                 call la_stgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                           dummy1,1,ifst,ilst,work(n*n*2 + 1),lwork - 2*n*n,ierr)
                 if (ierr > 0) then
                    ! ill-conditioned problem - swap rejected.
                    dif(ks) = zero
                 else
                    ! reordering successful, solve generalized sylvester
                    ! equation for r and l,
                               ! a22 * r - l * a11 = a12
                               ! b22 * r - l * b11 = b12,
                    ! and compute estimate of difl((a11,b11), (a22, b22)).
                    n1 = 1
                    if (work(2) /= zero) n1 = 2
                    n2 = n - n1
                    if (n2 == 0) then
                       dif(ks) = cond
                    else
                       i = n*n + 1
                       iz = 2*n*n + 1
                       call la_stgsyl('N',difdri,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),work(iz + 1),lwork - 2*n*n,iwork,ierr)
                       if (pair) dif(ks) = min(max(one,alprqt)*dif(ks),cond)
                    end if
                 end if
                 if (pair) dif(ks + 1) = dif(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_stgsna
     !> DTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B) in
     !> generalized real Schur canonical form (or of any matrix pair
     !> (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where
     !> Z**T denotes the transpose of Z.
     !> (A, B) must be in generalized real Schur form (as returned by DGGES),
     !> i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
     !> blocks. B is upper triangular.

     pure subroutine la_dtgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_dp,only:zero,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           real(dp),intent(out) :: dif(*),s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: difdri = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,iz,k,ks,lwmin,n1,n2
           real(dp) :: alphai,alphar,alprqt,beta,c1,c2,cond,eps,lnrm,rnrm,root1,root2, &
                     scale,smlnum,tmpii,tmpir,tmpri,tmprr,uhav,uhavi,uhbv,uhbvi
           ! Local Arrays
           real(dp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 pair = .false.
                 do k = 1,n
                    if (pair) then
                       pair = .false.
                    else
                       if (k < n) then
                          if (a(k + 1,k) == zero) then
                             if (select(k)) m = m + 1
                          else
                             pair = .true.
                             if (select(k) .or. select(k + 1)) m = m + 2
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*(n + 2) + 16
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              ! determine whether a(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_20
              else
                 if (k < n) pair = a(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_20
                 else
                    if (.not. select(k)) cycle loop_20
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (pair) then
                    ! complex eigenvalue pair.
                    rnrm = la_dlapy2(la_dnrm2(n,vr(1,ks),1),la_dnrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_dlapy2(la_dnrm2(n,vl(1,ks),1),la_dnrm2(n,vl( &
                              1,ks + 1),1))
                    call la_dgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    tmprr = la_ddot(n,work,1,vl(1,ks),1)
                    tmpri = la_ddot(n,work,1,vl(1,ks + 1),1)
                    call la_dgemv('N',n,n,one,a,lda,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_ddot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_ddot(n,work,1,vl(1,ks),1)
                    uhav = tmprr + tmpii
                    uhavi = tmpir - tmpri
                    call la_dgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    tmprr = la_ddot(n,work,1,vl(1,ks),1)
                    tmpri = la_ddot(n,work,1,vl(1,ks + 1),1)
                    call la_dgemv('N',n,n,one,b,ldb,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_ddot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_ddot(n,work,1,vl(1,ks),1)
                    uhbv = tmprr + tmpii
                    uhbvi = tmpir - tmpri
                    uhav = la_dlapy2(uhav,uhavi)
                    uhbv = la_dlapy2(uhbv,uhbvi)
                    cond = la_dlapy2(uhav,uhbv)
                    s(ks) = cond/(rnrm*lnrm)
                    s(ks + 1) = s(ks)
                 else
                    ! real eigenvalue.
                    rnrm = la_dnrm2(n,vr(1,ks),1)
                    lnrm = la_dnrm2(n,vl(1,ks),1)
                    call la_dgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    uhav = la_ddot(n,work,1,vl(1,ks),1)
                    call la_dgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    uhbv = la_ddot(n,work,1,vl(1,ks),1)
                    cond = la_dlapy2(uhav,uhbv)
                    if (cond == zero) then
                       s(ks) = -one
                    else
                       s(ks) = cond/(rnrm*lnrm)
                    end if
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_dlapy2(a(1,1),b(1,1))
                    cycle loop_20
                 end if
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvectors.
                 if (pair) then
                    ! copy the  2-by 2 pencil beginning at (a(k,k), b(k, k)).
                    ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_dlag2(work,2,work(5),2,smlnum*eps,beta,dummy1(1), &
                              alphar,dummy(1),alphai)
                    alprqt = one
                    c1 = two*(alphar*alphar + alphai*alphai + beta*beta)
                    c2 = four*beta*beta*alphai*alphai
                    root1 = c1 + sqrt(c1*c1 - 4.0_dp*c2)
                    root2 = c2/root1
                    root1 = root1/two
                    cond = min(sqrt(root1),sqrt(root2))
                 end if
                 ! copy the matrix (a, b) to the array work and swap the
                 ! diagonal block beginning at a(k,k) to the (1,1) position.
                 call la_dlacpy('FULL',n,n,a,lda,work,n)
                 call la_dlacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                 ifst = k
                 ilst = 1
                 call la_dtgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                           dummy1,1,ifst,ilst,work(n*n*2 + 1),lwork - 2*n*n,ierr)
                 if (ierr > 0) then
                    ! ill-conditioned problem - swap rejected.
                    dif(ks) = zero
                 else
                    ! reordering successful, solve generalized sylvester
                    ! equation for r and l,
                               ! a22 * r - l * a11 = a12
                               ! b22 * r - l * b11 = b12,
                    ! and compute estimate of difl((a11,b11), (a22, b22)).
                    n1 = 1
                    if (work(2) /= zero) n1 = 2
                    n2 = n - n1
                    if (n2 == 0) then
                       dif(ks) = cond
                    else
                       i = n*n + 1
                       iz = 2*n*n + 1
                       call la_dtgsyl('N',difdri,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),work(iz + 1),lwork - 2*n*n,iwork,ierr)
                       if (pair) dif(ks) = min(max(one,alprqt)*dif(ks),cond)
                    end if
                 end if
                 if (pair) dif(ks + 1) = dif(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_dtgsna
#ifdef LA_WITH_XDP
     !> XTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B) in
     !> generalized real Schur canonical form (or of any matrix pair
     !> (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where
     !> Z**T denotes the transpose of Z.
     !> (A, B) must be in generalized real Schur form (as returned by XGGES),
     !> i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
     !> blocks. B is upper triangular.

     pure subroutine la_xtgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_xdp,only:zero,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           real(xdp),intent(out) :: dif(*),s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: difdri = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,iz,k,ks,lwmin,n1,n2
           real(xdp) :: alphai,alphar,alprqt,beta,c1,c2,cond,eps,lnrm,rnrm,root1,root2, &
                     scale,smlnum,tmpii,tmpir,tmpri,tmprr,uhav,uhavi,uhbv,uhbvi
           ! Local Arrays
           real(xdp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 pair = .false.
                 do k = 1,n
                    if (pair) then
                       pair = .false.
                    else
                       if (k < n) then
                          if (a(k + 1,k) == zero) then
                             if (select(k)) m = m + 1
                          else
                             pair = .true.
                             if (select(k) .or. select(k + 1)) m = m + 2
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*(n + 2) + 16
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              ! determine whether a(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_20
              else
                 if (k < n) pair = a(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_20
                 else
                    if (.not. select(k)) cycle loop_20
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (pair) then
                    ! complex eigenvalue pair.
                    rnrm = la_xlapy2(la_xnrm2(n,vr(1,ks),1),la_xnrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_xlapy2(la_xnrm2(n,vl(1,ks),1),la_xnrm2(n,vl( &
                              1,ks + 1),1))
                    call la_xgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    tmprr = la_xdot(n,work,1,vl(1,ks),1)
                    tmpri = la_xdot(n,work,1,vl(1,ks + 1),1)
                    call la_xgemv('N',n,n,one,a,lda,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_xdot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_xdot(n,work,1,vl(1,ks),1)
                    uhav = tmprr + tmpii
                    uhavi = tmpir - tmpri
                    call la_xgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    tmprr = la_xdot(n,work,1,vl(1,ks),1)
                    tmpri = la_xdot(n,work,1,vl(1,ks + 1),1)
                    call la_xgemv('N',n,n,one,b,ldb,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_xdot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_xdot(n,work,1,vl(1,ks),1)
                    uhbv = tmprr + tmpii
                    uhbvi = tmpir - tmpri
                    uhav = la_xlapy2(uhav,uhavi)
                    uhbv = la_xlapy2(uhbv,uhbvi)
                    cond = la_xlapy2(uhav,uhbv)
                    s(ks) = cond/(rnrm*lnrm)
                    s(ks + 1) = s(ks)
                 else
                    ! real eigenvalue.
                    rnrm = la_xnrm2(n,vr(1,ks),1)
                    lnrm = la_xnrm2(n,vl(1,ks),1)
                    call la_xgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    uhav = la_xdot(n,work,1,vl(1,ks),1)
                    call la_xgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    uhbv = la_xdot(n,work,1,vl(1,ks),1)
                    cond = la_xlapy2(uhav,uhbv)
                    if (cond == zero) then
                       s(ks) = -one
                    else
                       s(ks) = cond/(rnrm*lnrm)
                    end if
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_xlapy2(a(1,1),b(1,1))
                    cycle loop_20
                 end if
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvectors.
                 if (pair) then
                    ! copy the  2-by 2 pencil beginning at (a(k,k), b(k, k)).
                    ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_xlag2(work,2,work(5),2,smlnum*eps,beta,dummy1(1), &
                              alphar,dummy(1),alphai)
                    alprqt = one
                    c1 = two*(alphar*alphar + alphai*alphai + beta*beta)
                    c2 = four*beta*beta*alphai*alphai
                    root1 = c1 + sqrt(c1*c1 - 4.0_xdp*c2)
                    root2 = c2/root1
                    root1 = root1/two
                    cond = min(sqrt(root1),sqrt(root2))
                 end if
                 ! copy the matrix (a, b) to the array work and swap the
                 ! diagonal block beginning at a(k,k) to the (1,1) position.
                 call la_xlacpy('FULL',n,n,a,lda,work,n)
                 call la_xlacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                 ifst = k
                 ilst = 1
                 call la_xtgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                           dummy1,1,ifst,ilst,work(n*n*2 + 1),lwork - 2*n*n,ierr)
                 if (ierr > 0) then
                    ! ill-conditioned problem - swap rejected.
                    dif(ks) = zero
                 else
                    ! reordering successful, solve generalized sylvester
                    ! equation for r and l,
                               ! a22 * r - l * a11 = a12
                               ! b22 * r - l * b11 = b12,
                    ! and compute estimate of difl((a11,b11), (a22, b22)).
                    n1 = 1
                    if (work(2) /= zero) n1 = 2
                    n2 = n - n1
                    if (n2 == 0) then
                       dif(ks) = cond
                    else
                       i = n*n + 1
                       iz = 2*n*n + 1
                       call la_xtgsyl('N',difdri,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),work(iz + 1),lwork - 2*n*n,iwork,ierr)
                       if (pair) dif(ks) = min(max(one,alprqt)*dif(ks),cond)
                    end if
                 end if
                 if (pair) dif(ks + 1) = dif(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_xtgsna
#endif
#ifdef LA_WITH_QP
     !> QTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B) in
     !> generalized real Schur canonical form (or of any matrix pair
     !> (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where
     !> Z**T denotes the transpose of Z.
     !> (A, B) must be in generalized real Schur form (as returned by QGGES),
     !> i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
     !> blocks. B is upper triangular.

     pure subroutine la_qtgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_qp,only:zero,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           real(qp),intent(out) :: dif(*),s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: difdri = 3

           ! Local Scalars
           logical(lk) :: lquery,pair,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,iz,k,ks,lwmin,n1,n2
           real(qp) :: alphai,alphar,alprqt,beta,c1,c2,cond,eps,lnrm,rnrm,root1,root2, &
                     scale,smlnum,tmpii,tmpir,tmpri,tmprr,uhav,uhavi,uhbv,uhbvi
           ! Local Arrays
           real(qp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 pair = .false.
                 do k = 1,n
                    if (pair) then
                       pair = .false.
                    else
                       if (k < n) then
                          if (a(k + 1,k) == zero) then
                             if (select(k)) m = m + 1
                          else
                             pair = .true.
                             if (select(k) .or. select(k + 1)) m = m + 2
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*(n + 2) + 16
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              ! determine whether a(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_20
              else
                 if (k < n) pair = a(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_20
                 else
                    if (.not. select(k)) cycle loop_20
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (pair) then
                    ! complex eigenvalue pair.
                    rnrm = la_qlapy2(la_qnrm2(n,vr(1,ks),1),la_qnrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_qlapy2(la_qnrm2(n,vl(1,ks),1),la_qnrm2(n,vl( &
                              1,ks + 1),1))
                    call la_qgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    tmprr = la_qdot(n,work,1,vl(1,ks),1)
                    tmpri = la_qdot(n,work,1,vl(1,ks + 1),1)
                    call la_qgemv('N',n,n,one,a,lda,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_qdot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_qdot(n,work,1,vl(1,ks),1)
                    uhav = tmprr + tmpii
                    uhavi = tmpir - tmpri
                    call la_qgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    tmprr = la_qdot(n,work,1,vl(1,ks),1)
                    tmpri = la_qdot(n,work,1,vl(1,ks + 1),1)
                    call la_qgemv('N',n,n,one,b,ldb,vr(1,ks + 1),1,zero,work,1)

                    tmpii = la_qdot(n,work,1,vl(1,ks + 1),1)
                    tmpir = la_qdot(n,work,1,vl(1,ks),1)
                    uhbv = tmprr + tmpii
                    uhbvi = tmpir - tmpri
                    uhav = la_qlapy2(uhav,uhavi)
                    uhbv = la_qlapy2(uhbv,uhbvi)
                    cond = la_qlapy2(uhav,uhbv)
                    s(ks) = cond/(rnrm*lnrm)
                    s(ks + 1) = s(ks)
                 else
                    ! real eigenvalue.
                    rnrm = la_qnrm2(n,vr(1,ks),1)
                    lnrm = la_qnrm2(n,vl(1,ks),1)
                    call la_qgemv('N',n,n,one,a,lda,vr(1,ks),1,zero,work,1)

                    uhav = la_qdot(n,work,1,vl(1,ks),1)
                    call la_qgemv('N',n,n,one,b,ldb,vr(1,ks),1,zero,work,1)

                    uhbv = la_qdot(n,work,1,vl(1,ks),1)
                    cond = la_qlapy2(uhav,uhbv)
                    if (cond == zero) then
                       s(ks) = -one
                    else
                       s(ks) = cond/(rnrm*lnrm)
                    end if
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_qlapy2(a(1,1),b(1,1))
                    cycle loop_20
                 end if
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvectors.
                 if (pair) then
                    ! copy the  2-by 2 pencil beginning at (a(k,k), b(k, k)).
                    ! compute the eigenvalue(s) at position k.
                    work(1) = a(k,k)
                    work(2) = a(k + 1,k)
                    work(3) = a(k,k + 1)
                    work(4) = a(k + 1,k + 1)
                    work(5) = b(k,k)
                    work(6) = b(k + 1,k)
                    work(7) = b(k,k + 1)
                    work(8) = b(k + 1,k + 1)
                    call la_qlag2(work,2,work(5),2,smlnum*eps,beta,dummy1(1), &
                              alphar,dummy(1),alphai)
                    alprqt = one
                    c1 = two*(alphar*alphar + alphai*alphai + beta*beta)
                    c2 = four*beta*beta*alphai*alphai
                    root1 = c1 + sqrt(c1*c1 - 4.0_qp*c2)
                    root2 = c2/root1
                    root1 = root1/two
                    cond = min(sqrt(root1),sqrt(root2))
                 end if
                 ! copy the matrix (a, b) to the array work and swap the
                 ! diagonal block beginning at a(k,k) to the (1,1) position.
                 call la_qlacpy('FULL',n,n,a,lda,work,n)
                 call la_qlacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                 ifst = k
                 ilst = 1
                 call la_qtgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                           dummy1,1,ifst,ilst,work(n*n*2 + 1),lwork - 2*n*n,ierr)
                 if (ierr > 0) then
                    ! ill-conditioned problem - swap rejected.
                    dif(ks) = zero
                 else
                    ! reordering successful, solve generalized sylvester
                    ! equation for r and l,
                               ! a22 * r - l * a11 = a12
                               ! b22 * r - l * b11 = b12,
                    ! and compute estimate of difl((a11,b11), (a22, b22)).
                    n1 = 1
                    if (work(2) /= zero) n1 = 2
                    n2 = n - n1
                    if (n2 == 0) then
                       dif(ks) = cond
                    else
                       i = n*n + 1
                       iz = 2*n*n + 1
                       call la_qtgsyl('N',difdri,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),work(iz + 1),lwork - 2*n*n,iwork,ierr)
                       if (pair) dif(ks) = min(max(one,alprqt)*dif(ks),cond)
                    end if
                 end if
                 if (pair) dif(ks + 1) = dif(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_qtgsna
#endif

     !> CTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of complex matrices (S,P), where S and P are upper triangular.
     !> Matrix pairs of this type are produced by the generalized Schur
     !> factorization of a complex matrix pair (A,B):
     !> A = Q*S*Z**H,  B = Q*P*Z**H
     !> as computed by CGGHRD + CHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal elements of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the unitary factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_ctgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,rwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: p(ldp,*),s(lds,*)
           complex(sp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: compl,compr,ilall,ilback,ilbbad,ilcomp,lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,im,iside,isrc,j,je,jr
           real(sp) :: acoefa,acoeff,anorm,ascale,bcoefa,big,bignum,bnorm,bscale,dmin, &
                     safmin,sbeta,scale,small,temp,ulp,xmax
           complex(sp) :: bcoeff,ca,cb,d,salpha,sum,suma,sumb,x
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,min,real
           ! Statement Functions
           real(sp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=sp)) + abs(aimag(x))
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors
           if (.not. ilall) then
              im = 0
              do j = 1,n
                 if (select(j)) im = im + 1
              end do
           else
              im = n
           end if
           ! check diagonal of b
           ilbbad = .false.
           do j = 1,n
              if (aimag(p(j,j)) /= zero) ilbbad = .true.
           end do
           if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('CTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_slamch('SAFE MINIMUM')
           big = one/safmin
           call la_slabad(safmin,big)
           ulp = la_slamch('EPSILON')*la_slamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part of a and b to check for possible overflow in the triangular
           ! solver.
           anorm = abs1(s(1,1))
           bnorm = abs1(p(1,1))
           rwork(1) = zero
           rwork(n + 1) = zero
           do j = 2,n
              rwork(j) = zero
              rwork(n + j) = zero
              do i = 1,j - 1
                 rwork(j) = rwork(j) + abs1(s(i,j))
                 rwork(n + j) = rwork(n + j) + abs1(p(i,j))
              end do
              anorm = max(anorm,rwork(j) + abs1(s(j,j)))
              bnorm = max(bnorm,rwork(n + j) + abs1(p(j,j)))
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              loop_140: do je = 1,n
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig + 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=sp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vl(jr,ieig) = czero
                       end do
                       vl(ieig,ieig) = cone
                       cycle loop_140
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                         ! h
                       ! y  ( a a - b b ) = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=sp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=sp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                    ! h
                    ! triangular solve of  (a a - b b)  y = 0
                                            ! h
                    ! (rowwise in  (a a - b b) , or columnwise in a a - b b)
                    loop_100: do j = je + 1,n
                       ! compute
                             ! j-1
                       ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                             ! k=je
                       ! (scale if necessary)
                       temp = one/xmax
                       if (acoefa*rwork(j) + bcoefa*rwork(n + j) > bignum*temp) then
                          do jr = je,j - 1
                             work(jr) = temp*work(jr)
                          end do
                          xmax = one
                       end if
                       suma = czero
                       sumb = czero
                       do jr = je,j - 1
                          suma = suma + conjg(s(jr,j))*work(jr)
                          sumb = sumb + conjg(p(jr,j))*work(jr)
                       end do
                       sum = acoeff*suma - conjg(bcoeff)*sumb
                       ! form x(j) = - sum / conjg( a*s(j,j) - b*p(j,j) )
                       ! with scaling and perturbation of the denominator
                       d = conjg(acoeff*s(j,j) - bcoeff*p(j,j))
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=sp)
                       if (abs1(d) < one) then
                          if (abs1(sum) >= bignum*abs1(d)) then
                             temp = one/abs1(sum)
                             do jr = je,j - 1
                                work(jr) = temp*work(jr)
                             end do
                             xmax = temp*xmax
                             sum = temp*sum
                          end if
                       end if
                       work(j) = la_cladiv(-sum,d)
                       xmax = max(xmax,abs1(work(j)))
                    end do loop_100
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_cgemv('N',n,n + 1 - je,cone,vl(1,je),ldvl,work(je),1, &
                                 czero,work(n + 1),1)
                       isrc = 2
                       ibeg = 1
                    else
                       isrc = 1
                       ibeg = je
                    end if
                    ! copy and scale eigenvector into column of vl
                    xmax = zero
                    do jr = ibeg,n
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = ibeg,n
                          vl(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       ibeg = n + 1
                    end if
                    do jr = 1,ibeg - 1
                       vl(jr,ieig) = czero
                    end do
                 end if
              end do loop_140
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              loop_250: do je = n,1,-1
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig - 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=sp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vr(jr,ieig) = czero
                       end do
                       vr(ieig,ieig) = cone
                       cycle loop_250
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                    ! ( a a - b b ) x  = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=sp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=sp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                    ! triangular solve of  (a a - b b) x = 0  (columnwise)
                    ! work(1:j-1) contains sums w,
                    ! work(j+1:je) contains x
                    do jr = 1,je - 1
                       work(jr) = acoeff*s(jr,je) - bcoeff*p(jr,je)
                    end do
                    work(je) = cone
                    loop_210: do j = je - 1,1,-1
                       ! form x(j) := - w(j) / d
                       ! with scaling and perturbation of the denominator
                       d = acoeff*s(j,j) - bcoeff*p(j,j)
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=sp)
                       if (abs1(d) < one) then
                          if (abs1(work(j)) >= bignum*abs1(d)) then
                             temp = one/abs1(work(j))
                             do jr = 1,je
                                work(jr) = temp*work(jr)
                             end do
                          end if
                       end if
                       work(j) = la_cladiv(-work(j),d)
                       if (j > 1) then
                          ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                          if (abs1(work(j)) > one) then
                             temp = one/abs1(work(j))
                             if (acoefa*rwork(j) + bcoefa*rwork(n + j) >= bignum*temp) then
                                do jr = 1,je
                                   work(jr) = temp*work(jr)
                                end do
                             end if
                          end if
                          ca = acoeff*work(j)
                          cb = bcoeff*work(j)
                          do jr = 1,j - 1
                             work(jr) = work(jr) + ca*s(jr,j) - cb*p(jr,j)
                          end do
                       end if
                    end do loop_210
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_cgemv('N',n,je,cone,vr,ldvr,work,1,czero,work(n + 1), &
                                 1)
                       isrc = 2
                       iend = n
                    else
                       isrc = 1
                       iend = je
                    end if
                    ! copy and scale eigenvector into column of vr
                    xmax = zero
                    do jr = 1,iend
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = 1,iend
                          vr(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       iend = 0
                    end if
                    do jr = iend + 1,n
                       vr(jr,ieig) = czero
                    end do
                 end if
              end do loop_250
           end if
           return
     end subroutine la_ctgevc
     !> ZTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of complex matrices (S,P), where S and P are upper triangular.
     !> Matrix pairs of this type are produced by the generalized Schur
     !> factorization of a complex matrix pair (A,B):
     !> A = Q*S*Z**H,  B = Q*P*Z**H
     !> as computed by ZGGHRD + ZHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal elements of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the unitary factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_ztgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,rwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: p(ldp,*),s(lds,*)
           complex(dp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: compl,compr,ilall,ilback,ilbbad,ilcomp,lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,im,iside,isrc,j,je,jr
           real(dp) :: acoefa,acoeff,anorm,ascale,bcoefa,big,bignum,bnorm,bscale,dmin, &
                     safmin,sbeta,scale,small,temp,ulp,xmax
           complex(dp) :: bcoeff,ca,cb,d,salpha,sum,suma,sumb,x
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(dp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=dp)) + abs(aimag(x))
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors
           if (.not. ilall) then
              im = 0
              do j = 1,n
                 if (select(j)) im = im + 1
              end do
           else
              im = n
           end if
           ! check diagonal of b
           ilbbad = .false.
           do j = 1,n
              if (aimag(p(j,j)) /= zero) ilbbad = .true.
           end do
           if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('ZTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_dlamch('SAFE MINIMUM')
           big = one/safmin
           call la_dlabad(safmin,big)
           ulp = la_dlamch('EPSILON')*la_dlamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part of a and b to check for possible overflow in the triangular
           ! solver.
           anorm = abs1(s(1,1))
           bnorm = abs1(p(1,1))
           rwork(1) = zero
           rwork(n + 1) = zero
           do j = 2,n
              rwork(j) = zero
              rwork(n + j) = zero
              do i = 1,j - 1
                 rwork(j) = rwork(j) + abs1(s(i,j))
                 rwork(n + j) = rwork(n + j) + abs1(p(i,j))
              end do
              anorm = max(anorm,rwork(j) + abs1(s(j,j)))
              bnorm = max(bnorm,rwork(n + j) + abs1(p(j,j)))
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              loop_140: do je = 1,n
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig + 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=dp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vl(jr,ieig) = czero
                       end do
                       vl(ieig,ieig) = cone
                       cycle loop_140
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                         ! h
                       ! y  ( a a - b b ) = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=dp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=dp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                    ! h
                    ! triangular solve of  (a a - b b)  y = 0
                                            ! h
                    ! (rowwise in  (a a - b b) , or columnwise in a a - b b)
                    loop_100: do j = je + 1,n
                       ! compute
                             ! j-1
                       ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                             ! k=je
                       ! (scale if necessary)
                       temp = one/xmax
                       if (acoefa*rwork(j) + bcoefa*rwork(n + j) > bignum*temp) then
                          do jr = je,j - 1
                             work(jr) = temp*work(jr)
                          end do
                          xmax = one
                       end if
                       suma = czero
                       sumb = czero
                       do jr = je,j - 1
                          suma = suma + conjg(s(jr,j))*work(jr)
                          sumb = sumb + conjg(p(jr,j))*work(jr)
                       end do
                       sum = acoeff*suma - conjg(bcoeff)*sumb
                       ! form x(j) = - sum / conjg( a*s(j,j) - b*p(j,j) )
                       ! with scaling and perturbation of the denominator
                       d = conjg(acoeff*s(j,j) - bcoeff*p(j,j))
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=dp)
                       if (abs1(d) < one) then
                          if (abs1(sum) >= bignum*abs1(d)) then
                             temp = one/abs1(sum)
                             do jr = je,j - 1
                                work(jr) = temp*work(jr)
                             end do
                             xmax = temp*xmax
                             sum = temp*sum
                          end if
                       end if
                       work(j) = la_zladiv(-sum,d)
                       xmax = max(xmax,abs1(work(j)))
                    end do loop_100
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_zgemv('N',n,n + 1 - je,cone,vl(1,je),ldvl,work(je),1, &
                                 czero,work(n + 1),1)
                       isrc = 2
                       ibeg = 1
                    else
                       isrc = 1
                       ibeg = je
                    end if
                    ! copy and scale eigenvector into column of vl
                    xmax = zero
                    do jr = ibeg,n
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = ibeg,n
                          vl(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       ibeg = n + 1
                    end if
                    do jr = 1,ibeg - 1
                       vl(jr,ieig) = czero
                    end do
                 end if
              end do loop_140
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              loop_250: do je = n,1,-1
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig - 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=dp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vr(jr,ieig) = czero
                       end do
                       vr(ieig,ieig) = cone
                       cycle loop_250
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                    ! ( a a - b b ) x  = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=dp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=dp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                    ! triangular solve of  (a a - b b) x = 0  (columnwise)
                    ! work(1:j-1) contains sums w,
                    ! work(j+1:je) contains x
                    do jr = 1,je - 1
                       work(jr) = acoeff*s(jr,je) - bcoeff*p(jr,je)
                    end do
                    work(je) = cone
                    loop_210: do j = je - 1,1,-1
                       ! form x(j) := - w(j) / d
                       ! with scaling and perturbation of the denominator
                       d = acoeff*s(j,j) - bcoeff*p(j,j)
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=dp)
                       if (abs1(d) < one) then
                          if (abs1(work(j)) >= bignum*abs1(d)) then
                             temp = one/abs1(work(j))
                             do jr = 1,je
                                work(jr) = temp*work(jr)
                             end do
                          end if
                       end if
                       work(j) = la_zladiv(-work(j),d)
                       if (j > 1) then
                          ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                          if (abs1(work(j)) > one) then
                             temp = one/abs1(work(j))
                             if (acoefa*rwork(j) + bcoefa*rwork(n + j) >= bignum*temp) then
                                do jr = 1,je
                                   work(jr) = temp*work(jr)
                                end do
                             end if
                          end if
                          ca = acoeff*work(j)
                          cb = bcoeff*work(j)
                          do jr = 1,j - 1
                             work(jr) = work(jr) + ca*s(jr,j) - cb*p(jr,j)
                          end do
                       end if
                    end do loop_210
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_zgemv('N',n,je,cone,vr,ldvr,work,1,czero,work(n + 1), &
                                 1)
                       isrc = 2
                       iend = n
                    else
                       isrc = 1
                       iend = je
                    end if
                    ! copy and scale eigenvector into column of vr
                    xmax = zero
                    do jr = 1,iend
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = 1,iend
                          vr(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       iend = 0
                    end if
                    do jr = iend + 1,n
                       vr(jr,ieig) = czero
                    end do
                 end if
              end do loop_250
           end if
           return
     end subroutine la_ztgevc
#ifdef LA_WITH_XDP
     !> YTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of complex matrices (S,P), where S and P are upper triangular.
     !> Matrix pairs of this type are produced by the generalized Schur
     !> factorization of a complex matrix pair (A,B):
     !> A = Q*S*Z**H,  B = Q*P*Z**H
     !> as computed by YGGHRD + YHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal elements of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the unitary factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_ytgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,rwork,info)
        use la_constants_xdp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(in) :: p(ldp,*),s(lds,*)
           complex(xdp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: compl,compr,ilall,ilback,ilbbad,ilcomp,lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,im,iside,isrc,j,je,jr
           real(xdp) :: acoefa,acoeff,anorm,ascale,bcoefa,big,bignum,bnorm,bscale,dmin, &
                     safmin,sbeta,scale,small,temp,ulp,xmax
           complex(xdp) :: bcoeff,ca,cb,d,salpha,sum,suma,sumb,x
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(xdp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=xdp)) + abs(aimag(x))
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('YTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors
           if (.not. ilall) then
              im = 0
              do j = 1,n
                 if (select(j)) im = im + 1
              end do
           else
              im = n
           end if
           ! check diagonal of b
           ilbbad = .false.
           do j = 1,n
              if (aimag(p(j,j)) /= zero) ilbbad = .true.
           end do
           if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('YTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_xlamch('SAFE MINIMUM')
           big = one/safmin
           call la_xlabad(safmin,big)
           ulp = la_xlamch('EPSILON')*la_xlamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part of a and b to check for possible overflow in the triangular
           ! solver.
           anorm = abs1(s(1,1))
           bnorm = abs1(p(1,1))
           rwork(1) = zero
           rwork(n + 1) = zero
           do j = 2,n
              rwork(j) = zero
              rwork(n + j) = zero
              do i = 1,j - 1
                 rwork(j) = rwork(j) + abs1(s(i,j))
                 rwork(n + j) = rwork(n + j) + abs1(p(i,j))
              end do
              anorm = max(anorm,rwork(j) + abs1(s(j,j)))
              bnorm = max(bnorm,rwork(n + j) + abs1(p(j,j)))
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              loop_140: do je = 1,n
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig + 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=xdp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vl(jr,ieig) = czero
                       end do
                       vl(ieig,ieig) = cone
                       cycle loop_140
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                         ! h
                       ! y  ( a a - b b ) = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=xdp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=xdp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                    ! h
                    ! triangular solve of  (a a - b b)  y = 0
                                            ! h
                    ! (rowwise in  (a a - b b) , or columnwise in a a - b b)
                    loop_100: do j = je + 1,n
                       ! compute
                             ! j-1
                       ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                             ! k=je
                       ! (scale if necessary)
                       temp = one/xmax
                       if (acoefa*rwork(j) + bcoefa*rwork(n + j) > bignum*temp) then
                          do jr = je,j - 1
                             work(jr) = temp*work(jr)
                          end do
                          xmax = one
                       end if
                       suma = czero
                       sumb = czero
                       do jr = je,j - 1
                          suma = suma + conjg(s(jr,j))*work(jr)
                          sumb = sumb + conjg(p(jr,j))*work(jr)
                       end do
                       sum = acoeff*suma - conjg(bcoeff)*sumb
                       ! form x(j) = - sum / conjg( a*s(j,j) - b*p(j,j) )
                       ! with scaling and perturbation of the denominator
                       d = conjg(acoeff*s(j,j) - bcoeff*p(j,j))
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=xdp)
                       if (abs1(d) < one) then
                          if (abs1(sum) >= bignum*abs1(d)) then
                             temp = one/abs1(sum)
                             do jr = je,j - 1
                                work(jr) = temp*work(jr)
                             end do
                             xmax = temp*xmax
                             sum = temp*sum
                          end if
                       end if
                       work(j) = la_yladiv(-sum,d)
                       xmax = max(xmax,abs1(work(j)))
                    end do loop_100
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_ygemv('N',n,n + 1 - je,cone,vl(1,je),ldvl,work(je),1, &
                                 czero,work(n + 1),1)
                       isrc = 2
                       ibeg = 1
                    else
                       isrc = 1
                       ibeg = je
                    end if
                    ! copy and scale eigenvector into column of vl
                    xmax = zero
                    do jr = ibeg,n
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = ibeg,n
                          vl(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       ibeg = n + 1
                    end if
                    do jr = 1,ibeg - 1
                       vl(jr,ieig) = czero
                    end do
                 end if
              end do loop_140
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              loop_250: do je = n,1,-1
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig - 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=xdp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vr(jr,ieig) = czero
                       end do
                       vr(ieig,ieig) = cone
                       cycle loop_250
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                    ! ( a a - b b ) x  = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=xdp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=xdp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                    ! triangular solve of  (a a - b b) x = 0  (columnwise)
                    ! work(1:j-1) contains sums w,
                    ! work(j+1:je) contains x
                    do jr = 1,je - 1
                       work(jr) = acoeff*s(jr,je) - bcoeff*p(jr,je)
                    end do
                    work(je) = cone
                    loop_210: do j = je - 1,1,-1
                       ! form x(j) := - w(j) / d
                       ! with scaling and perturbation of the denominator
                       d = acoeff*s(j,j) - bcoeff*p(j,j)
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=xdp)
                       if (abs1(d) < one) then
                          if (abs1(work(j)) >= bignum*abs1(d)) then
                             temp = one/abs1(work(j))
                             do jr = 1,je
                                work(jr) = temp*work(jr)
                             end do
                          end if
                       end if
                       work(j) = la_yladiv(-work(j),d)
                       if (j > 1) then
                          ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                          if (abs1(work(j)) > one) then
                             temp = one/abs1(work(j))
                             if (acoefa*rwork(j) + bcoefa*rwork(n + j) >= bignum*temp) then
                                do jr = 1,je
                                   work(jr) = temp*work(jr)
                                end do
                             end if
                          end if
                          ca = acoeff*work(j)
                          cb = bcoeff*work(j)
                          do jr = 1,j - 1
                             work(jr) = work(jr) + ca*s(jr,j) - cb*p(jr,j)
                          end do
                       end if
                    end do loop_210
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_ygemv('N',n,je,cone,vr,ldvr,work,1,czero,work(n + 1), &
                                 1)
                       isrc = 2
                       iend = n
                    else
                       isrc = 1
                       iend = je
                    end if
                    ! copy and scale eigenvector into column of vr
                    xmax = zero
                    do jr = 1,iend
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = 1,iend
                          vr(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       iend = 0
                    end if
                    do jr = iend + 1,n
                       vr(jr,ieig) = czero
                    end do
                 end if
              end do loop_250
           end if
           return
     end subroutine la_ytgevc
#endif
#ifdef LA_WITH_QP
     !> WTGEVC: computes some or all of the right and/or left eigenvectors of
     !> a pair of complex matrices (S,P), where S and P are upper triangular.
     !> Matrix pairs of this type are produced by the generalized Schur
     !> factorization of a complex matrix pair (A,B):
     !> A = Q*S*Z**H,  B = Q*P*Z**H
     !> as computed by WGGHRD + WHGEQZ.
     !> The right eigenvector x and the left eigenvector y of (S,P)
     !> corresponding to an eigenvalue w are defined by:
     !> S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
     !> where y**H denotes the conjugate tranpose of y.
     !> The eigenvalues are not input to this routine, but are computed
     !> directly from the diagonal elements of S and P.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
     !> where Z and Q are input matrices.
     !> If Q and Z are the unitary factors from the generalized Schur
     !> factorization of a matrix pair (A,B), then Z*X and Q*Y
     !> are the matrices of right and left eigenvectors of (A,B).

     pure subroutine la_wtgevc(side,howmny,select,n,s,lds,p,ldp,vl,ldvl,vr,ldvr, &
               mm,m,work,rwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldp,lds,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: p(ldp,*),s(lds,*)
           complex(qp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: compl,compr,ilall,ilback,ilbbad,ilcomp,lsa,lsb
           integer(ilp) :: i,ibeg,ieig,iend,ihwmny,im,iside,isrc,j,je,jr
           real(qp) :: acoefa,acoeff,anorm,ascale,bcoefa,big,bignum,bnorm,bscale,dmin, &
                     safmin,sbeta,scale,small,temp,ulp,xmax
           complex(qp) :: bcoeff,ca,cb,d,salpha,sum,suma,sumb,x
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(qp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=qp)) + abs(aimag(x))
           ! Executable Statements
           ! decode and test the input parameters
           if (la_lsame(howmny,'A')) then
              ihwmny = 1
              ilall = .true.
              ilback = .false.
           else if (la_lsame(howmny,'S')) then
              ihwmny = 2
              ilall = .false.
              ilback = .false.
           else if (la_lsame(howmny,'B')) then
              ihwmny = 3
              ilall = .true.
              ilback = .true.
           else
              ihwmny = -1
           end if
           if (la_lsame(side,'R')) then
              iside = 1
              compl = .false.
              compr = .true.
           else if (la_lsame(side,'L')) then
              iside = 2
              compl = .true.
              compr = .false.
           else if (la_lsame(side,'B')) then
              iside = 3
              compl = .true.
              compr = .true.
           else
              iside = -1
           end if
           info = 0
           if (iside < 0) then
              info = -1
           else if (ihwmny < 0) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lds < max(1,n)) then
              info = -6
           else if (ldp < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WTGEVC',-info)
              return
           end if
           ! count the number of eigenvectors
           if (.not. ilall) then
              im = 0
              do j = 1,n
                 if (select(j)) im = im + 1
              end do
           else
              im = n
           end if
           ! check diagonal of b
           ilbbad = .false.
           do j = 1,n
              if (aimag(p(j,j)) /= zero) ilbbad = .true.
           end do
           if (ilbbad) then
              info = -7
           else if (compl .and. ldvl < n .or. ldvl < 1) then
              info = -10
           else if (compr .and. ldvr < n .or. ldvr < 1) then
              info = -12
           else if (mm < im) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('WTGEVC',-info)
              return
           end if
           ! quick return if possible
           m = im
           if (n == 0) return
           ! machine constants
           safmin = la_qlamch('SAFE MINIMUM')
           big = one/safmin
           call la_qlabad(safmin,big)
           ulp = la_qlamch('EPSILON')*la_qlamch('BASE')
           small = safmin*n/ulp
           big = one/small
           bignum = one/(safmin*n)
           ! compute the 1-norm of each column of the strictly upper triangular
           ! part of a and b to check for possible overflow in the triangular
           ! solver.
           anorm = abs1(s(1,1))
           bnorm = abs1(p(1,1))
           rwork(1) = zero
           rwork(n + 1) = zero
           do j = 2,n
              rwork(j) = zero
              rwork(n + j) = zero
              do i = 1,j - 1
                 rwork(j) = rwork(j) + abs1(s(i,j))
                 rwork(n + j) = rwork(n + j) + abs1(p(i,j))
              end do
              anorm = max(anorm,rwork(j) + abs1(s(j,j)))
              bnorm = max(bnorm,rwork(n + j) + abs1(p(j,j)))
           end do
           ascale = one/max(anorm,safmin)
           bscale = one/max(bnorm,safmin)
           ! left eigenvectors
           if (compl) then
              ieig = 0
              ! main loop over eigenvalues
              loop_140: do je = 1,n
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig + 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=qp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vl(jr,ieig) = czero
                       end do
                       vl(ieig,ieig) = cone
                       cycle loop_140
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                         ! h
                       ! y  ( a a - b b ) = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=qp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=qp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                                                    ! h
                    ! triangular solve of  (a a - b b)  y = 0
                                            ! h
                    ! (rowwise in  (a a - b b) , or columnwise in a a - b b)
                    loop_100: do j = je + 1,n
                       ! compute
                             ! j-1
                       ! sum = sum  conjg( a*s(k,j) - b*p(k,j) )*x(k)
                             ! k=je
                       ! (scale if necessary)
                       temp = one/xmax
                       if (acoefa*rwork(j) + bcoefa*rwork(n + j) > bignum*temp) then
                          do jr = je,j - 1
                             work(jr) = temp*work(jr)
                          end do
                          xmax = one
                       end if
                       suma = czero
                       sumb = czero
                       do jr = je,j - 1
                          suma = suma + conjg(s(jr,j))*work(jr)
                          sumb = sumb + conjg(p(jr,j))*work(jr)
                       end do
                       sum = acoeff*suma - conjg(bcoeff)*sumb
                       ! form x(j) = - sum / conjg( a*s(j,j) - b*p(j,j) )
                       ! with scaling and perturbation of the denominator
                       d = conjg(acoeff*s(j,j) - bcoeff*p(j,j))
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=qp)
                       if (abs1(d) < one) then
                          if (abs1(sum) >= bignum*abs1(d)) then
                             temp = one/abs1(sum)
                             do jr = je,j - 1
                                work(jr) = temp*work(jr)
                             end do
                             xmax = temp*xmax
                             sum = temp*sum
                          end if
                       end if
                       work(j) = la_wladiv(-sum,d)
                       xmax = max(xmax,abs1(work(j)))
                    end do loop_100
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_wgemv('N',n,n + 1 - je,cone,vl(1,je),ldvl,work(je),1, &
                                 czero,work(n + 1),1)
                       isrc = 2
                       ibeg = 1
                    else
                       isrc = 1
                       ibeg = je
                    end if
                    ! copy and scale eigenvector into column of vl
                    xmax = zero
                    do jr = ibeg,n
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = ibeg,n
                          vl(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       ibeg = n + 1
                    end if
                    do jr = 1,ibeg - 1
                       vl(jr,ieig) = czero
                    end do
                 end if
              end do loop_140
           end if
           ! right eigenvectors
           if (compr) then
              ieig = im + 1
              ! main loop over eigenvalues
              loop_250: do je = n,1,-1
                 if (ilall) then
                    ilcomp = .true.
                 else
                    ilcomp = select(je)
                 end if
                 if (ilcomp) then
                    ieig = ieig - 1
                    if (abs1(s(je,je)) <= safmin .and. abs(real(p(je,je),KIND=qp)) &
                              <= safmin) then
                       ! singular matrix pencil -- return unit eigenvector
                       do jr = 1,n
                          vr(jr,ieig) = czero
                       end do
                       vr(ieig,ieig) = cone
                       cycle loop_250
                    end if
                    ! non-singular eigenvalue:
                    ! compute coefficients  a  and  b  in
                    ! ( a a - b b ) x  = 0
                    temp = one/max(abs1(s(je,je))*ascale,abs(real(p(je,je),KIND=qp)) &
                              *bscale,safmin)
                    salpha = (temp*s(je,je))*ascale
                    sbeta = (temp*real(p(je,je),KIND=qp))*bscale
                    acoeff = sbeta*ascale
                    bcoeff = salpha*bscale
                    ! scale to avoid underflow
                    lsa = abs(sbeta) >= safmin .and. abs(acoeff) < small
                    lsb = abs1(salpha) >= safmin .and. abs1(bcoeff) < small
                    scale = one
                    if (lsa) scale = (small/abs(sbeta))*min(anorm,big)
                    if (lsb) scale = max(scale, (small/abs1(salpha))*min(bnorm,big))

                    if (lsa .or. lsb) then
                       scale = min(scale,one/(safmin*max(one,abs(acoeff),abs1(bcoeff)) &
                                 ))
                       if (lsa) then
                          acoeff = ascale*(scale*sbeta)
                       else
                          acoeff = scale*acoeff
                       end if
                       if (lsb) then
                          bcoeff = bscale*(scale*salpha)
                       else
                          bcoeff = scale*bcoeff
                       end if
                    end if
                    acoefa = abs(acoeff)
                    bcoefa = abs1(bcoeff)
                    xmax = one
                    do jr = 1,n
                       work(jr) = czero
                    end do
                    work(je) = cone
                    dmin = max(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
                    ! triangular solve of  (a a - b b) x = 0  (columnwise)
                    ! work(1:j-1) contains sums w,
                    ! work(j+1:je) contains x
                    do jr = 1,je - 1
                       work(jr) = acoeff*s(jr,je) - bcoeff*p(jr,je)
                    end do
                    work(je) = cone
                    loop_210: do j = je - 1,1,-1
                       ! form x(j) := - w(j) / d
                       ! with scaling and perturbation of the denominator
                       d = acoeff*s(j,j) - bcoeff*p(j,j)
                       if (abs1(d) <= dmin) d = cmplx(dmin,KIND=qp)
                       if (abs1(d) < one) then
                          if (abs1(work(j)) >= bignum*abs1(d)) then
                             temp = one/abs1(work(j))
                             do jr = 1,je
                                work(jr) = temp*work(jr)
                             end do
                          end if
                       end if
                       work(j) = la_wladiv(-work(j),d)
                       if (j > 1) then
                          ! w = w + x(j)*(a s(*,j) - b p(*,j) ) with scaling
                          if (abs1(work(j)) > one) then
                             temp = one/abs1(work(j))
                             if (acoefa*rwork(j) + bcoefa*rwork(n + j) >= bignum*temp) then
                                do jr = 1,je
                                   work(jr) = temp*work(jr)
                                end do
                             end if
                          end if
                          ca = acoeff*work(j)
                          cb = bcoeff*work(j)
                          do jr = 1,j - 1
                             work(jr) = work(jr) + ca*s(jr,j) - cb*p(jr,j)
                          end do
                       end if
                    end do loop_210
                    ! back transform eigenvector if howmny='b'.
                    if (ilback) then
                       call la_wgemv('N',n,je,cone,vr,ldvr,work,1,czero,work(n + 1), &
                                 1)
                       isrc = 2
                       iend = n
                    else
                       isrc = 1
                       iend = je
                    end if
                    ! copy and scale eigenvector into column of vr
                    xmax = zero
                    do jr = 1,iend
                       xmax = max(xmax,abs1(work((isrc - 1)*n + jr)))
                    end do
                    if (xmax > safmin) then
                       temp = one/xmax
                       do jr = 1,iend
                          vr(jr,ieig) = temp*work((isrc - 1)*n + jr)
                       end do
                    else
                       iend = 0
                    end if
                    do jr = iend + 1,n
                       vr(jr,ieig) = czero
                    end do
                 end if
              end do loop_250
           end if
           return
     end subroutine la_wtgevc
#endif

     !> CTGEX2: swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
     !> in an upper triangular matrix pair (A, B) by an unitary equivalence
     !> transformation.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_ctgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,info)
        use la_constants_sp,only:czero,cone

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: twenty = 2.0e+1_sp
           integer(ilp),parameter :: ldst = 2
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,m
           real(sp) :: cq,cz,eps,sa,sb,scale,smlnum,sum,thresha,threshb
           complex(sp) :: cdum,f,g,sq,sz
           ! Local Arrays
           complex(sp) :: s(ldst,ldst),t(ldst,ldst),work(8)
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1) return
           m = ldst
           weak = .false.
           strong = .false.
           ! make a local copy of selected block in (a, b)
           call la_clacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_clacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute the threshold for testing the acceptance of swapping.
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           scale = real(czero,KIND=sp)
           sum = real(cone,KIND=sp)
           call la_clacpy('FULL',m,m,s,ldst,work,m)
           call la_clacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
           call la_classq(m*m,work,1,scale,sum)
           sa = scale*sqrt(sum)
           scale = real(czero,KIND=sp)
           sum = real(cone,KIND=sp)
           call la_classq(m*m,work(m*m + 1),1,scale,sum)
           sb = scale*sqrt(sum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*sa,smlnum)
           threshb = max(twenty*eps*sb,smlnum)
           ! compute unitary ql and rq that swap 1-by-1 and 1-by-1 blocks
           ! using givens rotations and perform the swap tentatively.
           f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
           g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
           sa = abs(s(2,2))*abs(t(1,1))
           sb = abs(s(1,1))*abs(t(2,2))
           call la_clartg(g,f,cz,sz,cdum)
           sz = -sz
           call la_crot(2,s(1,1),1,s(1,2),1,cz,conjg(sz))
           call la_crot(2,t(1,1),1,t(1,2),1,cz,conjg(sz))
           if (sa >= sb) then
              call la_clartg(s(1,1),s(2,1),cq,sq,cdum)
           else
              call la_clartg(t(1,1),t(2,1),cq,sq,cdum)
           end if
           call la_crot(2,s(1,1),ldst,s(2,1),ldst,cq,sq)
           call la_crot(2,t(1,1),ldst,t(2,1),ldst,cq,sq)
           ! weak stability test: |s21| <= o(eps f-norm((a)))
                                ! and  |t21| <= o(eps f-norm((b)))
           weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
           if (.not. weak) go to 20
           if (wands) then
              ! strong stability test:
                 ! f-norm((a-ql**h*s*qr, b-ql**h*t*qr)) <= o(eps*f-norm((a, b)))
              call la_clacpy('FULL',m,m,s,ldst,work,m)
              call la_clacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
              call la_crot(2,work,1,work(3),1,cz,-conjg(sz))
              call la_crot(2,work(5),1,work(7),1,cz,-conjg(sz))
              call la_crot(2,work,2,work(2),2,cq,-sq)
              call la_crot(2,work(5),2,work(6),2,cq,-sq)
              do i = 1,2
                 work(i) = work(i) - a(j1 + i - 1,j1)
                 work(i + 2) = work(i + 2) - a(j1 + i - 1,j1 + 1)
                 work(i + 4) = work(i + 4) - b(j1 + i - 1,j1)
                 work(i + 6) = work(i + 6) - b(j1 + i - 1,j1 + 1)
              end do
              scale = real(czero,KIND=sp)
              sum = real(cone,KIND=sp)
              call la_classq(m*m,work,1,scale,sum)
              sa = scale*sqrt(sum)
              scale = real(czero,KIND=sp)
              sum = real(cone,KIND=sp)
              call la_classq(m*m,work(m*m + 1),1,scale,sum)
              sb = scale*sqrt(sum)
              strong = sa <= thresha .and. sb <= threshb
              if (.not. strong) go to 20
           end if
           ! if the swap is accepted ("weakly" and "strongly"), apply the
           ! equivalence transformations to the original matrix pair (a,b)
           call la_crot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,cz,conjg(sz))
           call la_crot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,cz,conjg(sz))
           call la_crot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,cq,sq)
           call la_crot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,cq,sq)
           ! set  n1 by n2 (2,1) blocks to 0
           a(j1 + 1,j1) = czero
           b(j1 + 1,j1) = czero
           ! accumulate transformations into q and z if requested.
           if (wantz) call la_crot(n,z(1,j1),1,z(1,j1 + 1),1,cz,conjg(sz))

           if (wantq) call la_crot(n,q(1,j1),1,q(1,j1 + 1),1,cq,conjg(sq))

           ! exit with info = 0 if swap was successfully performed.
           return
           ! exit with info = 1 if swap was rejected.
           20 continue
           info = 1
           return
     end subroutine la_ctgex2
     !> ZTGEX2: swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
     !> in an upper triangular matrix pair (A, B) by an unitary equivalence
     !> transformation.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_ztgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,info)
        use la_constants_dp,only:czero,cone

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: twenty = 2.0e+1_dp
           integer(ilp),parameter :: ldst = 2
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,m
           real(dp) :: cq,cz,eps,sa,sb,scale,smlnum,sum,thresha,threshb
           complex(dp) :: cdum,f,g,sq,sz
           ! Local Arrays
           complex(dp) :: s(ldst,ldst),t(ldst,ldst),work(8)
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1) return
           m = ldst
           weak = .false.
           strong = .false.
           ! make a local copy of selected block in (a, b)
           call la_zlacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_zlacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute the threshold for testing the acceptance of swapping.
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           scale = real(czero,KIND=dp)
           sum = real(cone,KIND=dp)
           call la_zlacpy('FULL',m,m,s,ldst,work,m)
           call la_zlacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
           call la_zlassq(m*m,work,1,scale,sum)
           sa = scale*sqrt(sum)
           scale = real(czero,KIND=dp)
           sum = real(cone,KIND=dp)
           call la_zlassq(m*m,work(m*m + 1),1,scale,sum)
           sb = scale*sqrt(sum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*sa,smlnum)
           threshb = max(twenty*eps*sb,smlnum)
           ! compute unitary ql and rq that swap 1-by-1 and 1-by-1 blocks
           ! using givens rotations and perform the swap tentatively.
           f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
           g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
           sa = abs(s(2,2))*abs(t(1,1))
           sb = abs(s(1,1))*abs(t(2,2))
           call la_zlartg(g,f,cz,sz,cdum)
           sz = -sz
           call la_zrot(2,s(1,1),1,s(1,2),1,cz,conjg(sz))
           call la_zrot(2,t(1,1),1,t(1,2),1,cz,conjg(sz))
           if (sa >= sb) then
              call la_zlartg(s(1,1),s(2,1),cq,sq,cdum)
           else
              call la_zlartg(t(1,1),t(2,1),cq,sq,cdum)
           end if
           call la_zrot(2,s(1,1),ldst,s(2,1),ldst,cq,sq)
           call la_zrot(2,t(1,1),ldst,t(2,1),ldst,cq,sq)
           ! weak stability test: |s21| <= o(eps f-norm((a)))
                                ! and  |t21| <= o(eps f-norm((b)))
           weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
           if (.not. weak) go to 20
           if (wands) then
              ! strong stability test:
                 ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                 ! and
                 ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
              call la_zlacpy('FULL',m,m,s,ldst,work,m)
              call la_zlacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
              call la_zrot(2,work,1,work(3),1,cz,-conjg(sz))
              call la_zrot(2,work(5),1,work(7),1,cz,-conjg(sz))
              call la_zrot(2,work,2,work(2),2,cq,-sq)
              call la_zrot(2,work(5),2,work(6),2,cq,-sq)
              do i = 1,2
                 work(i) = work(i) - a(j1 + i - 1,j1)
                 work(i + 2) = work(i + 2) - a(j1 + i - 1,j1 + 1)
                 work(i + 4) = work(i + 4) - b(j1 + i - 1,j1)
                 work(i + 6) = work(i + 6) - b(j1 + i - 1,j1 + 1)
              end do
              scale = real(czero,KIND=dp)
              sum = real(cone,KIND=dp)
              call la_zlassq(m*m,work,1,scale,sum)
              sa = scale*sqrt(sum)
              scale = real(czero,KIND=dp)
              sum = real(cone,KIND=dp)
              call la_zlassq(m*m,work(m*m + 1),1,scale,sum)
              sb = scale*sqrt(sum)
              strong = sa <= thresha .and. sb <= threshb
              if (.not. strong) go to 20
           end if
           ! if the swap is accepted ("weakly" and "strongly"), apply the
           ! equivalence transformations to the original matrix pair (a,b)
           call la_zrot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,cz,conjg(sz))
           call la_zrot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,cz,conjg(sz))
           call la_zrot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,cq,sq)
           call la_zrot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,cq,sq)
           ! set  n1 by n2 (2,1) blocks to 0
           a(j1 + 1,j1) = czero
           b(j1 + 1,j1) = czero
           ! accumulate transformations into q and z if requested.
           if (wantz) call la_zrot(n,z(1,j1),1,z(1,j1 + 1),1,cz,conjg(sz))

           if (wantq) call la_zrot(n,q(1,j1),1,q(1,j1 + 1),1,cq,conjg(sq))

           ! exit with info = 0 if swap was successfully performed.
           return
           ! exit with info = 1 if swap was rejected.
           20 continue
           info = 1
           return
     end subroutine la_ztgex2
#ifdef LA_WITH_XDP
     !> YTGEX2: swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
     !> in an upper triangular matrix pair (A, B) by an unitary equivalence
     !> transformation.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_ytgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,info)
        use la_constants_xdp,only:czero,cone

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: twenty = 2.0e+1_xdp
           integer(ilp),parameter :: ldst = 2
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,m
           real(xdp) :: cq,cz,eps,sa,sb,scale,smlnum,sum,thresha,threshb
           complex(xdp) :: cdum,f,g,sq,sz
           ! Local Arrays
           complex(xdp) :: s(ldst,ldst),t(ldst,ldst),work(8)
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1) return
           m = ldst
           weak = .false.
           strong = .false.
           ! make a local copy of selected block in (a, b)
           call la_ylacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_ylacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute the threshold for testing the acceptance of swapping.
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           scale = real(czero,KIND=xdp)
           sum = real(cone,KIND=xdp)
           call la_ylacpy('FULL',m,m,s,ldst,work,m)
           call la_ylacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
           call la_ylassq(m*m,work,1,scale,sum)
           sa = scale*sqrt(sum)
           scale = real(czero,KIND=xdp)
           sum = real(cone,KIND=xdp)
           call la_ylassq(m*m,work(m*m + 1),1,scale,sum)
           sb = scale*sqrt(sum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*sa,smlnum)
           threshb = max(twenty*eps*sb,smlnum)
           ! compute unitary ql and rq that swap 1-by-1 and 1-by-1 blocks
           ! using givens rotations and perform the swap tentatively.
           f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
           g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
           sa = abs(s(2,2))*abs(t(1,1))
           sb = abs(s(1,1))*abs(t(2,2))
           call la_ylartg(g,f,cz,sz,cdum)
           sz = -sz
           call la_yrot(2,s(1,1),1,s(1,2),1,cz,conjg(sz))
           call la_yrot(2,t(1,1),1,t(1,2),1,cz,conjg(sz))
           if (sa >= sb) then
              call la_ylartg(s(1,1),s(2,1),cq,sq,cdum)
           else
              call la_ylartg(t(1,1),t(2,1),cq,sq,cdum)
           end if
           call la_yrot(2,s(1,1),ldst,s(2,1),ldst,cq,sq)
           call la_yrot(2,t(1,1),ldst,t(2,1),ldst,cq,sq)
           ! weak stability test: |s21| <= o(eps f-norm((a)))
                                ! and  |t21| <= o(eps f-norm((b)))
           weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
           if (.not. weak) go to 20
           if (wands) then
              ! strong stability test:
                 ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                 ! and
                 ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
              call la_ylacpy('FULL',m,m,s,ldst,work,m)
              call la_ylacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
              call la_yrot(2,work,1,work(3),1,cz,-conjg(sz))
              call la_yrot(2,work(5),1,work(7),1,cz,-conjg(sz))
              call la_yrot(2,work,2,work(2),2,cq,-sq)
              call la_yrot(2,work(5),2,work(6),2,cq,-sq)
              do i = 1,2
                 work(i) = work(i) - a(j1 + i - 1,j1)
                 work(i + 2) = work(i + 2) - a(j1 + i - 1,j1 + 1)
                 work(i + 4) = work(i + 4) - b(j1 + i - 1,j1)
                 work(i + 6) = work(i + 6) - b(j1 + i - 1,j1 + 1)
              end do
              scale = real(czero,KIND=xdp)
              sum = real(cone,KIND=xdp)
              call la_ylassq(m*m,work,1,scale,sum)
              sa = scale*sqrt(sum)
              scale = real(czero,KIND=xdp)
              sum = real(cone,KIND=xdp)
              call la_ylassq(m*m,work(m*m + 1),1,scale,sum)
              sb = scale*sqrt(sum)
              strong = sa <= thresha .and. sb <= threshb
              if (.not. strong) go to 20
           end if
           ! if the swap is accepted ("weakly" and "strongly"), apply the
           ! equivalence transformations to the original matrix pair (a,b)
           call la_yrot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,cz,conjg(sz))
           call la_yrot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,cz,conjg(sz))
           call la_yrot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,cq,sq)
           call la_yrot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,cq,sq)
           ! set  n1 by n2 (2,1) blocks to 0
           a(j1 + 1,j1) = czero
           b(j1 + 1,j1) = czero
           ! accumulate transformations into q and z if requested.
           if (wantz) call la_yrot(n,z(1,j1),1,z(1,j1 + 1),1,cz,conjg(sz))

           if (wantq) call la_yrot(n,q(1,j1),1,q(1,j1 + 1),1,cq,conjg(sq))

           ! exit with info = 0 if swap was successfully performed.
           return
           ! exit with info = 1 if swap was rejected.
           20 continue
           info = 1
           return
     end subroutine la_ytgex2
#endif
#ifdef LA_WITH_QP
     !> WTGEX2: swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
     !> in an upper triangular matrix pair (A, B) by an unitary equivalence
     !> transformation.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_wtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,j1,info)
        use la_constants_qp,only:czero,cone

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,lda,ldb,ldq,ldz,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: twenty = 2.0e+1_qp
           integer(ilp),parameter :: ldst = 2
           logical(lk),parameter :: wands = .true.

           ! Local Scalars
           logical(lk) :: strong,weak
           integer(ilp) :: i,m
           real(qp) :: cq,cz,eps,sa,sb,scale,smlnum,sum,thresha,threshb
           complex(qp) :: cdum,f,g,sq,sz
           ! Local Arrays
           complex(qp) :: s(ldst,ldst),t(ldst,ldst),work(8)
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,max,sqrt
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 1) return
           m = ldst
           weak = .false.
           strong = .false.
           ! make a local copy of selected block in (a, b)
           call la_wlacpy('FULL',m,m,a(j1,j1),lda,s,ldst)
           call la_wlacpy('FULL',m,m,b(j1,j1),ldb,t,ldst)
           ! compute the threshold for testing the acceptance of swapping.
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           scale = real(czero,KIND=qp)
           sum = real(cone,KIND=qp)
           call la_wlacpy('FULL',m,m,s,ldst,work,m)
           call la_wlacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
           call la_wlassq(m*m,work,1,scale,sum)
           sa = scale*sqrt(sum)
           scale = real(czero,KIND=qp)
           sum = real(cone,KIND=qp)
           call la_wlassq(m*m,work(m*m + 1),1,scale,sum)
           sb = scale*sqrt(sum)
           ! thres has been changed from
              ! thresh = max( ten*eps*sa, smlnum )
           ! to
              ! thresh = max( twenty*eps*sa, smlnum )
           ! on 04/01/10.
           ! "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
           ! jim demmel and guillaume revy. see forum post 1783.
           thresha = max(twenty*eps*sa,smlnum)
           threshb = max(twenty*eps*sb,smlnum)
           ! compute unitary ql and rq that swap 1-by-1 and 1-by-1 blocks
           ! using givens rotations and perform the swap tentatively.
           f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
           g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
           sa = abs(s(2,2))*abs(t(1,1))
           sb = abs(s(1,1))*abs(t(2,2))
           call la_wlartg(g,f,cz,sz,cdum)
           sz = -sz
           call la_wrot(2,s(1,1),1,s(1,2),1,cz,conjg(sz))
           call la_wrot(2,t(1,1),1,t(1,2),1,cz,conjg(sz))
           if (sa >= sb) then
              call la_wlartg(s(1,1),s(2,1),cq,sq,cdum)
           else
              call la_wlartg(t(1,1),t(2,1),cq,sq,cdum)
           end if
           call la_wrot(2,s(1,1),ldst,s(2,1),ldst,cq,sq)
           call la_wrot(2,t(1,1),ldst,t(2,1),ldst,cq,sq)
           ! weak stability test: |s21| <= o(eps f-norm((a)))
                                ! and  |t21| <= o(eps f-norm((b)))
           weak = abs(s(2,1)) <= thresha .and. abs(t(2,1)) <= threshb
           if (.not. weak) go to 20
           if (wands) then
              ! strong stability test:
                 ! f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
                 ! and
                 ! f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
              call la_wlacpy('FULL',m,m,s,ldst,work,m)
              call la_wlacpy('FULL',m,m,t,ldst,work(m*m + 1),m)
              call la_wrot(2,work,1,work(3),1,cz,-conjg(sz))
              call la_wrot(2,work(5),1,work(7),1,cz,-conjg(sz))
              call la_wrot(2,work,2,work(2),2,cq,-sq)
              call la_wrot(2,work(5),2,work(6),2,cq,-sq)
              do i = 1,2
                 work(i) = work(i) - a(j1 + i - 1,j1)
                 work(i + 2) = work(i + 2) - a(j1 + i - 1,j1 + 1)
                 work(i + 4) = work(i + 4) - b(j1 + i - 1,j1)
                 work(i + 6) = work(i + 6) - b(j1 + i - 1,j1 + 1)
              end do
              scale = real(czero,KIND=qp)
              sum = real(cone,KIND=qp)
              call la_wlassq(m*m,work,1,scale,sum)
              sa = scale*sqrt(sum)
              scale = real(czero,KIND=qp)
              sum = real(cone,KIND=qp)
              call la_wlassq(m*m,work(m*m + 1),1,scale,sum)
              sb = scale*sqrt(sum)
              strong = sa <= thresha .and. sb <= threshb
              if (.not. strong) go to 20
           end if
           ! if the swap is accepted ("weakly" and "strongly"), apply the
           ! equivalence transformations to the original matrix pair (a,b)
           call la_wrot(j1 + 1,a(1,j1),1,a(1,j1 + 1),1,cz,conjg(sz))
           call la_wrot(j1 + 1,b(1,j1),1,b(1,j1 + 1),1,cz,conjg(sz))
           call la_wrot(n - j1 + 1,a(j1,j1),lda,a(j1 + 1,j1),lda,cq,sq)
           call la_wrot(n - j1 + 1,b(j1,j1),ldb,b(j1 + 1,j1),ldb,cq,sq)
           ! set  n1 by n2 (2,1) blocks to 0
           a(j1 + 1,j1) = czero
           b(j1 + 1,j1) = czero
           ! accumulate transformations into q and z if requested.
           if (wantz) call la_wrot(n,z(1,j1),1,z(1,j1 + 1),1,cz,conjg(sz))

           if (wantq) call la_wrot(n,q(1,j1),1,q(1,j1 + 1),1,cq,conjg(sq))

           ! exit with info = 0 if swap was successfully performed.
           return
           ! exit with info = 1 if swap was rejected.
           20 continue
           info = 1
           return
     end subroutine la_wtgex2
#endif

     !> CTGEXC: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A,B), using an unitary equivalence transformation
     !> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
     !> row index IFST is moved to row ILST.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_ctgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ifst,lda,ldb,ldq,ldz,n
           integer(ilp),intent(inout) :: ilst
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: here
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('CTGEXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           if (ifst == ilst) return
           if (ifst < ilst) then
              here = ifst
              10 continue
              ! swap with next one below
              call la_ctgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here + 1
              if (here < ilst) go to 10
              here = here - 1
           else
              here = ifst - 1
              20 continue
              ! swap with next one above
              call la_ctgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here - 1
              if (here >= ilst) go to 20
              here = here + 1
           end if
           ilst = here
           return
     end subroutine la_ctgexc
     !> ZTGEXC: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A,B), using an unitary equivalence transformation
     !> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
     !> row index IFST is moved to row ILST.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_ztgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ifst,lda,ldb,ldq,ldz,n
           integer(ilp),intent(inout) :: ilst
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: here
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('ZTGEXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           if (ifst == ilst) return
           if (ifst < ilst) then
              here = ifst
              10 continue
              ! swap with next one below
              call la_ztgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here + 1
              if (here < ilst) go to 10
              here = here - 1
           else
              here = ifst - 1
              20 continue
              ! swap with next one above
              call la_ztgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here - 1
              if (here >= ilst) go to 20
              here = here + 1
           end if
           ilst = here
           return
     end subroutine la_ztgexc
#ifdef LA_WITH_XDP
     !> YTGEXC: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A,B), using an unitary equivalence transformation
     !> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
     !> row index IFST is moved to row ILST.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_ytgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ifst,lda,ldb,ldq,ldz,n
           integer(ilp),intent(inout) :: ilst
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: here
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('YTGEXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           if (ifst == ilst) return
           if (ifst < ilst) then
              here = ifst
              10 continue
              ! swap with next one below
              call la_ytgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here + 1
              if (here < ilst) go to 10
              here = here - 1
           else
              here = ifst - 1
              20 continue
              ! swap with next one above
              call la_ytgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here - 1
              if (here >= ilst) go to 20
              here = here + 1
           end if
           ilst = here
           return
     end subroutine la_ytgexc
#endif
#ifdef LA_WITH_QP
     !> WTGEXC: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A,B), using an unitary equivalence transformation
     !> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
     !> row index IFST is moved to row ILST.
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.
     !> Optionally, the matrices Q and Z of generalized Schur vectors are
     !> updated.
     !> Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
     !> Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H

     pure subroutine la_wtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,ifst,ilst, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ifst,lda,ldb,ldq,ldz,n
           integer(ilp),intent(inout) :: ilst
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: here
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test input arguments.
           info = 0
           if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldq < 1 .or. wantq .and. (ldq < max(1,n))) then
              info = -9
           else if (ldz < 1 .or. wantz .and. (ldz < max(1,n))) then
              info = -11
           else if (ifst < 1 .or. ifst > n) then
              info = -12
           else if (ilst < 1 .or. ilst > n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('WTGEXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           if (ifst == ilst) return
           if (ifst < ilst) then
              here = ifst
              10 continue
              ! swap with next one below
              call la_wtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here + 1
              if (here < ilst) go to 10
              here = here - 1
           else
              here = ifst - 1
              20 continue
              ! swap with next one above
              call la_wtgex2(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,here,info)

              if (info /= 0) then
                 ilst = here
                 return
              end if
              here = here - 1
              if (here >= ilst) go to 20
              here = here + 1
           end if
           ilst = here
           return
     end subroutine la_wtgexc
#endif

     !> CTGSY2: solves the generalized Sylvester equation
     !> A * R - L * B = scale *  C               (1)
     !> D * R - L * E = scale * F
     !> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
     !> (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Zx = scale * b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Ik is the identity matrix of size k and X**H is the transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R  + D**H * L   = scale * C           (3)
     !> R  * B**H + L  * E**H  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> = sigma_min(Z) using reverse communication with CLACON.
     !> CTGSY2 also (IJOB >= 1) contributes to the computation in CTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of two matrix pairs in
     !> CTGSYL.

     pure subroutine la_ctgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info
           real(sp),intent(inout) :: rdscal,rdsum
           real(sp),intent(out) :: scale
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(sp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldz = 2

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ierr,j,k
           real(sp) :: scaloc
           complex(sp) :: alpha
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           complex(sp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: cmplx,conjg,max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CTGSY2',-info)
              return
           end if
           if (notran) then
              ! solve (i, j) - system
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = m, m - 1, ..., 1; j = 1, 2, ..., n
              scale = one
              scaloc = one
              loop_30: do j = 1,n
                 loop_20: do i = m,1,-1
                    ! build 2 by 2 system
                    z(1,1) = a(i,i)
                    z(2,1) = d(i,i)
                    z(1,2) = -b(j,j)
                    z(2,2) = -e(j,j)
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z * x = rhs
                    call la_cgetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    if (ijob == 0) then
                       call la_cgesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_cscal(m,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                             call la_cscal(m,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                    else
                       call la_clatdf(ijob,ldz,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (i > 1) then
                       alpha = -rhs(1)
                       call la_caxpy(i - 1,alpha,a(1,i),1,c(1,j),1)
                       call la_caxpy(i - 1,alpha,d(1,i),1,f(1,j),1)
                    end if
                    if (j < n) then
                       call la_caxpy(n - j,rhs(2),b(j,j + 1),ldb,c(i,j + 1),ldc)

                       call la_caxpy(n - j,rhs(2),e(j,j + 1),lde,f(i,j + 1),ldf)

                    end if
                 end do loop_20
              end do loop_30
           else
              ! solve transposed (i, j) - system:
                 ! a(i, i)**h * r(i, j) + d(i, i)**h * l(j, j) = c(i, j)
                 ! r(i, i) * b(j, j) + l(i, j) * e(j, j)   = -f(i, j)
              ! for i = 1, 2, ..., m, j = n, n - 1, ..., 1
              scale = one
              scaloc = one
              loop_80: do i = 1,m
                 loop_70: do j = n,1,-1
                    ! build 2 by 2 system z**h
                    z(1,1) = conjg(a(i,i))
                    z(2,1) = -conjg(b(j,j))
                    z(1,2) = conjg(d(i,i))
                    z(2,2) = -conjg(e(j,j))
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z**h * x = rhs
                    call la_cgetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    call la_cgesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                    if (scaloc /= one) then
                       do k = 1,n
                          call la_cscal(m,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                          call la_cscal(m,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    do k = 1,j - 1
                       f(i,k) = f(i,k) + rhs(1)*conjg(b(k,j)) + rhs(2)*conjg(e(k, &
                                 j))
                    end do
                    do k = i + 1,m
                       c(k,j) = c(k,j) - conjg(a(i,k))*rhs(1) - conjg(d(i,k)) &
                                 *rhs(2)
                    end do
                 end do loop_70
              end do loop_80
           end if
           return
     end subroutine la_ctgsy2
     !> ZTGSY2: solves the generalized Sylvester equation
     !> A * R - L * B = scale * C               (1)
     !> D * R - L * E = scale * F
     !> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
     !> (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Zx = scale * b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Ik is the identity matrix of size k and X**H is the conjuguate transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R  + D**H * L   = scale * C           (3)
     !> R  * B**H + L  * E**H  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> = sigma_min(Z) using reverse communication with ZLACON.
     !> ZTGSY2 also (IJOB >= 1) contributes to the computation in ZTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of two matrix pairs in
     !> ZTGSYL.

     pure subroutine la_ztgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info
           real(dp),intent(inout) :: rdscal,rdsum
           real(dp),intent(out) :: scale
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(dp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldz = 2

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ierr,j,k
           real(dp) :: scaloc
           complex(dp) :: alpha
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           complex(dp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: cmplx,conjg,max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZTGSY2',-info)
              return
           end if
           if (notran) then
              ! solve (i, j) - system
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = m, m - 1, ..., 1; j = 1, 2, ..., n
              scale = one
              scaloc = one
              loop_30: do j = 1,n
                 loop_20: do i = m,1,-1
                    ! build 2 by 2 system
                    z(1,1) = a(i,i)
                    z(2,1) = d(i,i)
                    z(1,2) = -b(j,j)
                    z(2,2) = -e(j,j)
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z * x = rhs
                    call la_zgetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    if (ijob == 0) then
                       call la_zgesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_zscal(m,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                             call la_zscal(m,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                    else
                       call la_zlatdf(ijob,ldz,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (i > 1) then
                       alpha = -rhs(1)
                       call la_zaxpy(i - 1,alpha,a(1,i),1,c(1,j),1)
                       call la_zaxpy(i - 1,alpha,d(1,i),1,f(1,j),1)
                    end if
                    if (j < n) then
                       call la_zaxpy(n - j,rhs(2),b(j,j + 1),ldb,c(i,j + 1),ldc)

                       call la_zaxpy(n - j,rhs(2),e(j,j + 1),lde,f(i,j + 1),ldf)

                    end if
                 end do loop_20
              end do loop_30
           else
              ! solve transposed (i, j) - system:
                 ! a(i, i)**h * r(i, j) + d(i, i)**h * l(j, j) = c(i, j)
                 ! r(i, i) * b(j, j) + l(i, j) * e(j, j)   = -f(i, j)
              ! for i = 1, 2, ..., m, j = n, n - 1, ..., 1
              scale = one
              scaloc = one
              loop_80: do i = 1,m
                 loop_70: do j = n,1,-1
                    ! build 2 by 2 system z**h
                    z(1,1) = conjg(a(i,i))
                    z(2,1) = -conjg(b(j,j))
                    z(1,2) = conjg(d(i,i))
                    z(2,2) = -conjg(e(j,j))
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z**h * x = rhs
                    call la_zgetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    call la_zgesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                    if (scaloc /= one) then
                       do k = 1,n
                          call la_zscal(m,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                          call la_zscal(m,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    do k = 1,j - 1
                       f(i,k) = f(i,k) + rhs(1)*conjg(b(k,j)) + rhs(2)*conjg(e(k, &
                                 j))
                    end do
                    do k = i + 1,m
                       c(k,j) = c(k,j) - conjg(a(i,k))*rhs(1) - conjg(d(i,k)) &
                                 *rhs(2)
                    end do
                 end do loop_70
              end do loop_80
           end if
           return
     end subroutine la_ztgsy2
#ifdef LA_WITH_XDP
     !> YTGSY2: solves the generalized Sylvester equation
     !> A * R - L * B = scale * C               (1)
     !> D * R - L * E = scale * F
     !> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
     !> (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Zx = scale * b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Ik is the identity matrix of size k and X**H is the conjuguate transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R  + D**H * L   = scale * C           (3)
     !> R  * B**H + L  * E**H  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> = sigma_min(Z) using reverse communication with YLACON.
     !> YTGSY2 also (IJOB >= 1) contributes to the computation in YTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of two matrix pairs in
     !> YTGSYL.

     pure subroutine la_ytgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info
           real(xdp),intent(inout) :: rdscal,rdsum
           real(xdp),intent(out) :: scale
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(xdp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldz = 2

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ierr,j,k
           real(xdp) :: scaloc
           complex(xdp) :: alpha
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           complex(xdp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: cmplx,conjg,max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YTGSY2',-info)
              return
           end if
           if (notran) then
              ! solve (i, j) - system
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = m, m - 1, ..., 1; j = 1, 2, ..., n
              scale = one
              scaloc = one
              loop_30: do j = 1,n
                 loop_20: do i = m,1,-1
                    ! build 2 by 2 system
                    z(1,1) = a(i,i)
                    z(2,1) = d(i,i)
                    z(1,2) = -b(j,j)
                    z(2,2) = -e(j,j)
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z * x = rhs
                    call la_ygetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    if (ijob == 0) then
                       call la_ygesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                             call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                    else
                       call la_ylatdf(ijob,ldz,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (i > 1) then
                       alpha = -rhs(1)
                       call la_yaxpy(i - 1,alpha,a(1,i),1,c(1,j),1)
                       call la_yaxpy(i - 1,alpha,d(1,i),1,f(1,j),1)
                    end if
                    if (j < n) then
                       call la_yaxpy(n - j,rhs(2),b(j,j + 1),ldb,c(i,j + 1),ldc)

                       call la_yaxpy(n - j,rhs(2),e(j,j + 1),lde,f(i,j + 1),ldf)

                    end if
                 end do loop_20
              end do loop_30
           else
              ! solve transposed (i, j) - system:
                 ! a(i, i)**h * r(i, j) + d(i, i)**h * l(j, j) = c(i, j)
                 ! r(i, i) * b(j, j) + l(i, j) * e(j, j)   = -f(i, j)
              ! for i = 1, 2, ..., m, j = n, n - 1, ..., 1
              scale = one
              scaloc = one
              loop_80: do i = 1,m
                 loop_70: do j = n,1,-1
                    ! build 2 by 2 system z**h
                    z(1,1) = conjg(a(i,i))
                    z(2,1) = -conjg(b(j,j))
                    z(1,2) = conjg(d(i,i))
                    z(2,2) = -conjg(e(j,j))
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z**h * x = rhs
                    call la_ygetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    call la_ygesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                    if (scaloc /= one) then
                       do k = 1,n
                          call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                          call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    do k = 1,j - 1
                       f(i,k) = f(i,k) + rhs(1)*conjg(b(k,j)) + rhs(2)*conjg(e(k, &
                                 j))
                    end do
                    do k = i + 1,m
                       c(k,j) = c(k,j) - conjg(a(i,k))*rhs(1) - conjg(d(i,k)) &
                                 *rhs(2)
                    end do
                 end do loop_70
              end do loop_80
           end if
           return
     end subroutine la_ytgsy2
#endif
#ifdef LA_WITH_QP
     !> WTGSY2: solves the generalized Sylvester equation
     !> A * R - L * B = scale * C               (1)
     !> D * R - L * E = scale * F
     !> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
     !> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
     !> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
     !> (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
     !> scaling factor chosen to avoid overflow.
     !> In matrix notation solving equation (1) corresponds to solve
     !> Zx = scale * b, where Z is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Ik is the identity matrix of size k and X**H is the conjuguate transpose of X.
     !> kron(X, Y) is the Kronecker product between the matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R  + D**H * L   = scale * C           (3)
     !> R  * B**H + L  * E**H  = scale * -F
     !> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
     !> = sigma_min(Z) using reverse communication with WLACON.
     !> WTGSY2 also (IJOB >= 1) contributes to the computation in WTGSYL
     !> of an upper bound on the separation between to matrix pairs. Then
     !> the input (A, D), (B, E) are sub-pencils of two matrix pairs in
     !> WTGSYL.

     pure subroutine la_wtgsy2(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,rdsum,rdscal,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,m,n
           integer(ilp),intent(out) :: info
           real(qp),intent(inout) :: rdscal,rdsum
           real(qp),intent(out) :: scale
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(qp),intent(inout) :: c(ldc,*),f(ldf,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldz = 2

           ! Local Scalars
           logical(lk) :: notran
           integer(ilp) :: i,ierr,j,k
           real(qp) :: scaloc
           complex(qp) :: alpha
           ! Local Arrays
           integer(ilp) :: ipiv(ldz),jpiv(ldz)
           complex(qp) :: rhs(ldz),z(ldz,ldz)
           ! Intrinsic Functions
           intrinsic :: cmplx,conjg,max
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           ierr = 0
           notran = la_lsame(trans,'N')
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 2)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WTGSY2',-info)
              return
           end if
           if (notran) then
              ! solve (i, j) - system
                 ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                 ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
              ! for i = m, m - 1, ..., 1; j = 1, 2, ..., n
              scale = one
              scaloc = one
              loop_30: do j = 1,n
                 loop_20: do i = m,1,-1
                    ! build 2 by 2 system
                    z(1,1) = a(i,i)
                    z(2,1) = d(i,i)
                    z(1,2) = -b(j,j)
                    z(2,2) = -e(j,j)
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z * x = rhs
                    call la_wgetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    if (ijob == 0) then
                       call la_wgesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                       if (scaloc /= one) then
                          do k = 1,n
                             call la_wscal(m,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                             call la_wscal(m,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                    else
                       call la_wlatdf(ijob,ldz,z,ldz,rhs,rdsum,rdscal,ipiv,jpiv)

                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    if (i > 1) then
                       alpha = -rhs(1)
                       call la_waxpy(i - 1,alpha,a(1,i),1,c(1,j),1)
                       call la_waxpy(i - 1,alpha,d(1,i),1,f(1,j),1)
                    end if
                    if (j < n) then
                       call la_waxpy(n - j,rhs(2),b(j,j + 1),ldb,c(i,j + 1),ldc)

                       call la_waxpy(n - j,rhs(2),e(j,j + 1),lde,f(i,j + 1),ldf)

                    end if
                 end do loop_20
              end do loop_30
           else
              ! solve transposed (i, j) - system:
                 ! a(i, i)**h * r(i, j) + d(i, i)**h * l(j, j) = c(i, j)
                 ! r(i, i) * b(j, j) + l(i, j) * e(j, j)   = -f(i, j)
              ! for i = 1, 2, ..., m, j = n, n - 1, ..., 1
              scale = one
              scaloc = one
              loop_80: do i = 1,m
                 loop_70: do j = n,1,-1
                    ! build 2 by 2 system z**h
                    z(1,1) = conjg(a(i,i))
                    z(2,1) = -conjg(b(j,j))
                    z(1,2) = conjg(d(i,i))
                    z(2,2) = -conjg(e(j,j))
                    ! set up right hand side(s)
                    rhs(1) = c(i,j)
                    rhs(2) = f(i,j)
                    ! solve z**h * x = rhs
                    call la_wgetc2(ldz,z,ldz,ipiv,jpiv,ierr)
                    if (ierr > 0) info = ierr
                    call la_wgesc2(ldz,z,ldz,rhs,ipiv,jpiv,scaloc)
                    if (scaloc /= one) then
                       do k = 1,n
                          call la_wscal(m,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                          call la_wscal(m,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! unpack solution vector(s)
                    c(i,j) = rhs(1)
                    f(i,j) = rhs(2)
                    ! substitute r(i, j) and l(i, j) into remaining equation.
                    do k = 1,j - 1
                       f(i,k) = f(i,k) + rhs(1)*conjg(b(k,j)) + rhs(2)*conjg(e(k, &
                                 j))
                    end do
                    do k = i + 1,m
                       c(k,j) = c(k,j) - conjg(a(i,k))*rhs(1) - conjg(d(i,k)) &
                                 *rhs(2)
                    end do
                 end do loop_70
              end do loop_80
           end if
           return
     end subroutine la_wtgsy2
#endif

     !> CTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C            (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with complex entries. A, B, D and E are upper
     !> triangular (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1
     !> is an output scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z
     !> is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Here Ix is the identity matrix of size x and X**H is the conjugate
     !> transpose of X. Kron(X, Y) is the Kronecker product between the
     !> matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R + D**H * L = scale * C           (3)
     !> R * B**H + L * E**H = scale * -F
     !> This case (TRANS = 'C') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using CLACON.
     !> If IJOB >= 1, CTGSYL computes a Frobenius norm-based estimate of
     !> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z.
     !> This is a level-3 BLAS algorithm.

     pure subroutine la_ctgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(sp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           complex(sp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(sp),intent(inout) :: c(ldc,*),f(ldf,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_ccopy by calls to la_claset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,pq,q
           real(sp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: cmplx,max,real,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine  optimal block sizes mb and nb
           mb = la_ilaenv(2,'CTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'CTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_claset('F',m,n,czero,czero,c,ldc)
                 call la_claset('F',m,n,czero,czero,f,ldf)
              else if (ijob >= 1 .and. notran) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              ! use unblocked level 2 solver
              loop_30: do iround = 1,isolve
                 scale = one
                 dscale = zero
                 dsum = one
                 pq = m*n
                 call la_ctgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=sp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=sp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_clacpy('F',m,n,c,ldc,work,m)
                    call la_clacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_claset('F',m,n,czero,czero,c,ldc)
                    call la_claset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_clacpy('F',m,n,work,m,c,ldc)
                    call la_clacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j) - subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
                 pq = 0
                 scale = one
                 dscale = zero
                 dsum = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       call la_ctgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + mb*nb
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_cscal(m,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                             call la_cscal(m,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_cscal(is - 1,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                             call la_cscal(is - 1,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_cscal(m - ie,cmplx(scaloc,zero,KIND=sp),c(ie + 1,k), &
                                       1)
                             call la_cscal(m - ie,cmplx(scaloc,zero,KIND=sp),f(ie + 1,k), &
                                       1)
                          end do
                          do k = je + 1,n
                             call la_cscal(m,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                             call la_cscal(m,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i,j) and l(i,j) into remaining equation.
                       if (i > 1) then
                          call la_cgemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=sp),a( &
                           1,is),lda,c(is,js),ldc,cmplx(one,zero,KIND=sp),c(1,js), &
                                     ldc)
                          call la_cgemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=sp),d( &
                           1,is),ldd,c(is,js),ldc,cmplx(one,zero,KIND=sp),f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_cgemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=sp),f( &
                          is,js),ldf,b(js,je + 1),ldb,cmplx(one,zero,KIND=sp),c(is,je + 1 &
                                    ),ldc)
                          call la_cgemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=sp),f( &
                          is,js),ldf,e(js,je + 1),lde,cmplx(one,zero,KIND=sp),f(is,je + 1 &
                                    ),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=sp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=sp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_clacpy('F',m,n,c,ldc,work,m)
                    call la_clacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_claset('F',m,n,czero,czero,c,ldc)
                    call la_claset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_clacpy('F',m,n,work,m,c,ldc)
                    call la_clacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                  ! a(i, i)**h * r(i, j) + d(i, i)**h * l(i, j) = c(i, j)
                  ! r(i, j) * b(j, j)  + l(i, j) * e(j, j) = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_ctgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_cscal(m,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                          call la_cscal(m,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_cscal(is - 1,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                          call la_cscal(is - 1,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_cscal(m - ie,cmplx(scaloc,zero,KIND=sp),c(ie + 1,k),1)

                          call la_cscal(m - ie,cmplx(scaloc,zero,KIND=sp),f(ie + 1,k),1)

                       end do
                       do k = je + 1,n
                          call la_cscal(m,cmplx(scaloc,zero,KIND=sp),c(1,k),1)

                          call la_cscal(m,cmplx(scaloc,zero,KIND=sp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i,j) and l(i,j) into remaining equation.
                    if (j > p + 2) then
                       call la_cgemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=sp),c(is, &
                        js),ldc,b(1,js),ldb,cmplx(one,zero,KIND=sp),f(is,1),ldf)

                       call la_cgemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=sp),f(is, &
                        js),ldf,e(1,js),lde,cmplx(one,zero,KIND=sp),f(is,1),ldf)

                    end if
                    if (i < p) then
                       call la_cgemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=sp),a( &
                       is,ie + 1),lda,c(is,js),ldc,cmplx(one,zero,KIND=sp),c(ie + 1,js), &
                                 ldc)
                       call la_cgemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=sp),d( &
                       is,ie + 1),ldd,f(is,js),ldf,cmplx(one,zero,KIND=sp),c(ie + 1,js), &
                                 ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_ctgsyl
     !> ZTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C            (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with complex entries. A, B, D and E are upper
     !> triangular (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1
     !> is an output scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z
     !> is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Here Ix is the identity matrix of size x and X**H is the conjugate
     !> transpose of X. Kron(X, Y) is the Kronecker product between the
     !> matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R + D**H * L = scale * C           (3)
     !> R * B**H + L * E**H = scale * -F
     !> This case (TRANS = 'C') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using ZLACON.
     !> If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of
     !> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z.
     !> This is a level-3 BLAS algorithm.

     pure subroutine la_ztgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(dp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           complex(dp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(dp),intent(inout) :: c(ldc,*),f(ldf,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_ccopy by calls to la_claset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,pq,q
           real(dp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine  optimal block sizes mb and nb
           mb = la_ilaenv(2,'ZTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'ZTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_zlaset('F',m,n,czero,czero,c,ldc)
                 call la_zlaset('F',m,n,czero,czero,f,ldf)
              else if (ijob >= 1 .and. notran) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              ! use unblocked level 2 solver
              loop_30: do iround = 1,isolve
                 scale = one
                 dscale = zero
                 dsum = one
                 pq = m*n
                 call la_ztgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=dp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=dp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_zlacpy('F',m,n,c,ldc,work,m)
                    call la_zlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_zlaset('F',m,n,czero,czero,c,ldc)
                    call la_zlaset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_zlacpy('F',m,n,work,m,c,ldc)
                    call la_zlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j) - subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
                 pq = 0
                 scale = one
                 dscale = zero
                 dsum = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       call la_ztgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + mb*nb
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_zscal(m,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                             call la_zscal(m,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_zscal(is - 1,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                             call la_zscal(is - 1,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_zscal(m - ie,cmplx(scaloc,zero,KIND=dp),c(ie + 1,k), &
                                       1)
                             call la_zscal(m - ie,cmplx(scaloc,zero,KIND=dp),f(ie + 1,k), &
                                       1)
                          end do
                          do k = je + 1,n
                             call la_zscal(m,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                             call la_zscal(m,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i,j) and l(i,j) into remaining equation.
                       if (i > 1) then
                          call la_zgemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=dp),a( &
                           1,is),lda,c(is,js),ldc,cmplx(one,zero,KIND=dp),c(1,js), &
                                     ldc)
                          call la_zgemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=dp),d( &
                           1,is),ldd,c(is,js),ldc,cmplx(one,zero,KIND=dp),f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_zgemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=dp),f( &
                          is,js),ldf,b(js,je + 1),ldb,cmplx(one,zero,KIND=dp),c(is,je + 1 &
                                    ),ldc)
                          call la_zgemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=dp),f( &
                          is,js),ldf,e(js,je + 1),lde,cmplx(one,zero,KIND=dp),f(is,je + 1 &
                                    ),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=dp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=dp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_zlacpy('F',m,n,c,ldc,work,m)
                    call la_zlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_zlaset('F',m,n,czero,czero,c,ldc)
                    call la_zlaset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_zlacpy('F',m,n,work,m,c,ldc)
                    call la_zlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                  ! a(i, i)**h * r(i, j) + d(i, i)**h * l(i, j) = c(i, j)
                  ! r(i, j) * b(j, j)  + l(i, j) * e(j, j) = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_ztgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_zscal(m,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                          call la_zscal(m,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_zscal(is - 1,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                          call la_zscal(is - 1,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_zscal(m - ie,cmplx(scaloc,zero,KIND=dp),c(ie + 1,k),1)

                          call la_zscal(m - ie,cmplx(scaloc,zero,KIND=dp),f(ie + 1,k),1)

                       end do
                       do k = je + 1,n
                          call la_zscal(m,cmplx(scaloc,zero,KIND=dp),c(1,k),1)

                          call la_zscal(m,cmplx(scaloc,zero,KIND=dp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i,j) and l(i,j) into remaining equation.
                    if (j > p + 2) then
                       call la_zgemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=dp),c(is, &
                        js),ldc,b(1,js),ldb,cmplx(one,zero,KIND=dp),f(is,1),ldf)

                       call la_zgemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=dp),f(is, &
                        js),ldf,e(1,js),lde,cmplx(one,zero,KIND=dp),f(is,1),ldf)

                    end if
                    if (i < p) then
                       call la_zgemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=dp),a( &
                       is,ie + 1),lda,c(is,js),ldc,cmplx(one,zero,KIND=dp),c(ie + 1,js), &
                                 ldc)
                       call la_zgemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=dp),d( &
                       is,ie + 1),ldd,f(is,js),ldf,cmplx(one,zero,KIND=dp),c(ie + 1,js), &
                                 ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_ztgsyl
#ifdef LA_WITH_XDP
     !> YTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C            (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with complex entries. A, B, D and E are upper
     !> triangular (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1
     !> is an output scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z
     !> is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Here Ix is the identity matrix of size x and X**H is the conjugate
     !> transpose of X. Kron(X, Y) is the Kronecker product between the
     !> matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R + D**H * L = scale * C           (3)
     !> R * B**H + L * E**H = scale * -F
     !> This case (TRANS = 'C') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using YLACON.
     !> If IJOB >= 1, YTGSYL computes a Frobenius norm-based estimate of
     !> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z.
     !> This is a level-3 BLAS algorithm.

     pure subroutine la_ytgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(xdp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(xdp),intent(inout) :: c(ldc,*),f(ldf,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_zcopy by calls to la_zlaset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,pq,q
           real(xdp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine  optimal block sizes mb and nb
           mb = la_ilaenv(2,'YTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'YTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_ylaset('F',m,n,czero,czero,c,ldc)
                 call la_ylaset('F',m,n,czero,czero,f,ldf)
              else if (ijob >= 1 .and. notran) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              ! use unblocked level 2 solver
              loop_30: do iround = 1,isolve
                 scale = one
                 dscale = zero
                 dsum = one
                 pq = m*n
                 call la_ytgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=xdp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=xdp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_ylacpy('F',m,n,c,ldc,work,m)
                    call la_ylacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_ylaset('F',m,n,czero,czero,c,ldc)
                    call la_ylaset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_ylacpy('F',m,n,work,m,c,ldc)
                    call la_ylacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j) - subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
                 pq = 0
                 scale = one
                 dscale = zero
                 dsum = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       call la_ytgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + mb*nb
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                             call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_yscal(is - 1,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                             call la_yscal(is - 1,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_yscal(m - ie,cmplx(scaloc,zero,KIND=xdp),c(ie + 1,k), &
                                       1)
                             call la_yscal(m - ie,cmplx(scaloc,zero,KIND=xdp),f(ie + 1,k), &
                                       1)
                          end do
                          do k = je + 1,n
                             call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                             call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i,j) and l(i,j) into remaining equation.
                       if (i > 1) then
                          call la_ygemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=xdp),a( &
                           1,is),lda,c(is,js),ldc,cmplx(one,zero,KIND=xdp),c(1,js), &
                                     ldc)
                          call la_ygemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=xdp),d( &
                           1,is),ldd,c(is,js),ldc,cmplx(one,zero,KIND=xdp),f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_ygemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=xdp),f( &
                          is,js),ldf,b(js,je + 1),ldb,cmplx(one,zero,KIND=xdp),c(is,je + 1 &
                                    ),ldc)
                          call la_ygemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=xdp),f( &
                          is,js),ldf,e(js,je + 1),lde,cmplx(one,zero,KIND=xdp),f(is,je + 1 &
                                    ),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=xdp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=xdp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_ylacpy('F',m,n,c,ldc,work,m)
                    call la_ylacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_ylaset('F',m,n,czero,czero,c,ldc)
                    call la_ylaset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_ylacpy('F',m,n,work,m,c,ldc)
                    call la_ylacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                  ! a(i, i)**h * r(i, j) + d(i, i)**h * l(i, j) = c(i, j)
                  ! r(i, j) * b(j, j)  + l(i, j) * e(j, j) = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_ytgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                          call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_yscal(is - 1,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                          call la_yscal(is - 1,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_yscal(m - ie,cmplx(scaloc,zero,KIND=xdp),c(ie + 1,k),1)

                          call la_yscal(m - ie,cmplx(scaloc,zero,KIND=xdp),f(ie + 1,k),1)

                       end do
                       do k = je + 1,n
                          call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),c(1,k),1)

                          call la_yscal(m,cmplx(scaloc,zero,KIND=xdp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i,j) and l(i,j) into remaining equation.
                    if (j > p + 2) then
                       call la_ygemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=xdp),c(is, &
                        js),ldc,b(1,js),ldb,cmplx(one,zero,KIND=xdp),f(is,1),ldf)

                       call la_ygemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=xdp),f(is, &
                        js),ldf,e(1,js),lde,cmplx(one,zero,KIND=xdp),f(is,1),ldf)

                    end if
                    if (i < p) then
                       call la_ygemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=xdp),a( &
                       is,ie + 1),lda,c(is,js),ldc,cmplx(one,zero,KIND=xdp),c(ie + 1,js), &
                                 ldc)
                       call la_ygemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=xdp),d( &
                       is,ie + 1),ldd,f(is,js),ldf,cmplx(one,zero,KIND=xdp),c(ie + 1,js), &
                                 ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_ytgsyl
#endif
#ifdef LA_WITH_QP
     !> WTGSYL: solves the generalized Sylvester equation:
     !> A * R - L * B = scale * C            (1)
     !> D * R - L * E = scale * F
     !> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
     !> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
     !> respectively, with complex entries. A, B, D and E are upper
     !> triangular (i.e., (A,D) and (B,E) in generalized Schur form).
     !> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1
     !> is an output scaling factor chosen to avoid overflow.
     !> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z
     !> is defined as
     !> Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)
     !> [ kron(In, D)  -kron(E**H, Im) ],
     !> Here Ix is the identity matrix of size x and X**H is the conjugate
     !> transpose of X. Kron(X, Y) is the Kronecker product between the
     !> matrices X and Y.
     !> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b
     !> is solved for, which is equivalent to solve for R and L in
     !> A**H * R + D**H * L = scale * C           (3)
     !> R * B**H + L * E**H = scale * -F
     !> This case (TRANS = 'C') is used to compute an one-norm-based estimate
     !> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
     !> and (B,E), using WLACON.
     !> If IJOB >= 1, WTGSYL computes a Frobenius norm-based estimate of
     !> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
     !> reciprocal of the smallest singular value of Z.
     !> This is a level-3 BLAS algorithm.

     pure subroutine la_wtgsyl(trans,ijob,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
               ldf,scale,dif,work,lwork,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ijob,lda,ldb,ldc,ldd,lde,ldf,lwork,m,n
           integer(ilp),intent(out) :: info
           real(qp),intent(out) :: dif,scale
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           complex(qp),intent(in) :: a(lda,*),b(ldb,*),d(ldd,*),e(lde,*)
           complex(qp),intent(inout) :: c(ldc,*),f(ldf,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
        ! replaced various illegal calls to la_zcopy by calls to la_zlaset.
        ! sven hammarling, 1/5/02.

           ! Local Scalars
           logical(lk) :: lquery,notran
           integer(ilp) :: i,ie,ifunc,iround,is,isolve,j,je,js,k,linfo,lwmin,mb,nb, &
                     p,pq,q
           real(qp) :: dscale,dsum,scale2,scaloc
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max,sqrt
           ! Executable Statements
           ! decode and test input parameters
           info = 0
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -1
           else if (notran) then
              if ((ijob < 0) .or. (ijob > 4)) then
                 info = -2
              end if
           end if
           if (info == 0) then
              if (m <= 0) then
                 info = -3
              else if (n <= 0) then
                 info = -4
              else if (lda < max(1,m)) then
                 info = -6
              else if (ldb < max(1,n)) then
                 info = -8
              else if (ldc < max(1,m)) then
                 info = -10
              else if (ldd < max(1,m)) then
                 info = -12
              else if (lde < max(1,n)) then
                 info = -14
              else if (ldf < max(1,m)) then
                 info = -16
              end if
           end if
           if (info == 0) then
              if (notran) then
                 if (ijob == 1 .or. ijob == 2) then
                    lwmin = max(1,2*m*n)
                 else
                    lwmin = 1
                 end if
              else
                 lwmin = 1
              end if
              work(1) = lwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WTGSYL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              scale = 1
              if (notran) then
                 if (ijob /= 0) then
                    dif = 0
                 end if
              end if
              return
           end if
           ! determine  optimal block sizes mb and nb
           mb = la_ilaenv(2,'WTGSYL',trans,m,n,-1,-1)
           nb = la_ilaenv(5,'WTGSYL',trans,m,n,-1,-1)
           isolve = 1
           ifunc = 0
           if (notran) then
              if (ijob >= 3) then
                 ifunc = ijob - 2
                 call la_wlaset('F',m,n,czero,czero,c,ldc)
                 call la_wlaset('F',m,n,czero,czero,f,ldf)
              else if (ijob >= 1 .and. notran) then
                 isolve = 2
              end if
           end if
           if ((mb <= 1 .and. nb <= 1) .or. (mb >= m .and. nb >= n)) then
              ! use unblocked level 2 solver
              loop_30: do iround = 1,isolve
                 scale = one
                 dscale = zero
                 dsum = one
                 pq = m*n
                 call la_wtgsy2(trans,ifunc,m,n,a,lda,b,ldb,c,ldc,d,ldd,e,lde,f, &
                            ldf,scale,dsum,dscale,info)
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=qp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=qp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_wlacpy('F',m,n,c,ldc,work,m)
                    call la_wlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_wlaset('F',m,n,czero,czero,c,ldc)
                    call la_wlaset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_wlacpy('F',m,n,work,m,c,ldc)
                    call la_wlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_30
              return
           end if
           ! determine block structure of a
           p = 0
           i = 1
           40 continue
           if (i > m) go to 50
           p = p + 1
           iwork(p) = i
           i = i + mb
           if (i >= m) go to 50
           go to 40
           50 continue
           iwork(p + 1) = m + 1
           if (iwork(p) == iwork(p + 1)) p = p - 1
           ! determine block structure of b
           q = p + 1
           j = 1
           60 continue
           if (j > n) go to 70
           q = q + 1
           iwork(q) = j
           j = j + nb
           if (j >= n) go to 70
           go to 60
           70 continue
           iwork(q + 1) = n + 1
           if (iwork(q) == iwork(q + 1)) q = q - 1
           if (notran) then
              loop_150: do iround = 1,isolve
                 ! solve (i, j) - subsystem
                     ! a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
                     ! d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
                 ! for i = p, p - 1, ..., 1; j = 1, 2, ..., q
                 pq = 0
                 scale = one
                 dscale = zero
                 dsum = one
                 loop_130: do j = p + 2,q
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    loop_120: do i = p,1,-1
                       is = iwork(i)
                       ie = iwork(i + 1) - 1
                       mb = ie - is + 1
                       call la_wtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js), &
                       ldb,c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf, &
                                 scaloc,dsum,dscale,linfo)
                       if (linfo > 0) info = linfo
                       pq = pq + mb*nb
                       if (scaloc /= one) then
                          do k = 1,js - 1
                             call la_wscal(m,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                             call la_wscal(m,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_wscal(is - 1,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                             call la_wscal(is - 1,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                          end do
                          do k = js,je
                             call la_wscal(m - ie,cmplx(scaloc,zero,KIND=qp),c(ie + 1,k), &
                                       1)
                             call la_wscal(m - ie,cmplx(scaloc,zero,KIND=qp),f(ie + 1,k), &
                                       1)
                          end do
                          do k = je + 1,n
                             call la_wscal(m,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                             call la_wscal(m,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                          end do
                          scale = scale*scaloc
                       end if
                       ! substitute r(i,j) and l(i,j) into remaining equation.
                       if (i > 1) then
                          call la_wgemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=qp),a( &
                           1,is),lda,c(is,js),ldc,cmplx(one,zero,KIND=qp),c(1,js), &
                                     ldc)
                          call la_wgemm('N','N',is - 1,nb,mb,cmplx(-one,zero,KIND=qp),d( &
                           1,is),ldd,c(is,js),ldc,cmplx(one,zero,KIND=qp),f(1,js), &
                                     ldf)
                       end if
                       if (j < q) then
                          call la_wgemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=qp),f( &
                          is,js),ldf,b(js,je + 1),ldb,cmplx(one,zero,KIND=qp),c(is,je + 1 &
                                    ),ldc)
                          call la_wgemm('N','N',mb,n - je,nb,cmplx(one,zero,KIND=qp),f( &
                          is,js),ldf,e(js,je + 1),lde,cmplx(one,zero,KIND=qp),f(is,je + 1 &
                                    ),ldf)
                       end if
                    end do loop_120
                 end do loop_130
                 if (dscale /= zero) then
                    if (ijob == 1 .or. ijob == 3) then
                       dif = sqrt(real(2*m*n,KIND=qp))/(dscale*sqrt(dsum))
                    else
                       dif = sqrt(real(pq,KIND=qp))/(dscale*sqrt(dsum))
                    end if
                 end if
                 if (isolve == 2 .and. iround == 1) then
                    if (notran) then
                       ifunc = ijob
                    end if
                    scale2 = scale
                    call la_wlacpy('F',m,n,c,ldc,work,m)
                    call la_wlacpy('F',m,n,f,ldf,work(m*n + 1),m)
                    call la_wlaset('F',m,n,czero,czero,c,ldc)
                    call la_wlaset('F',m,n,czero,czero,f,ldf)
                 else if (isolve == 2 .and. iround == 2) then
                    call la_wlacpy('F',m,n,work,m,c,ldc)
                    call la_wlacpy('F',m,n,work(m*n + 1),m,f,ldf)
                    scale = scale2
                 end if
              end do loop_150
           else
              ! solve transposed (i, j)-subsystem
                  ! a(i, i)**h * r(i, j) + d(i, i)**h * l(i, j) = c(i, j)
                  ! r(i, j) * b(j, j)  + l(i, j) * e(j, j) = -f(i, j)
              ! for i = 1,2,..., p; j = q, q-1,..., 1
              scale = one
              loop_210: do i = 1,p
                 is = iwork(i)
                 ie = iwork(i + 1) - 1
                 mb = ie - is + 1
                 loop_200: do j = q,p + 2,-1
                    js = iwork(j)
                    je = iwork(j + 1) - 1
                    nb = je - js + 1
                    call la_wtgsy2(trans,ifunc,mb,nb,a(is,is),lda,b(js,js),ldb, &
                    c(is,js),ldc,d(is,is),ldd,e(js,js),lde,f(is,js),ldf,scaloc, &
                              dsum,dscale,linfo)
                    if (linfo > 0) info = linfo
                    if (scaloc /= one) then
                       do k = 1,js - 1
                          call la_wscal(m,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                          call la_wscal(m,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_wscal(is - 1,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                          call la_wscal(is - 1,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                       end do
                       do k = js,je
                          call la_wscal(m - ie,cmplx(scaloc,zero,KIND=qp),c(ie + 1,k),1)

                          call la_wscal(m - ie,cmplx(scaloc,zero,KIND=qp),f(ie + 1,k),1)

                       end do
                       do k = je + 1,n
                          call la_wscal(m,cmplx(scaloc,zero,KIND=qp),c(1,k),1)

                          call la_wscal(m,cmplx(scaloc,zero,KIND=qp),f(1,k),1)

                       end do
                       scale = scale*scaloc
                    end if
                    ! substitute r(i,j) and l(i,j) into remaining equation.
                    if (j > p + 2) then
                       call la_wgemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=qp),c(is, &
                        js),ldc,b(1,js),ldb,cmplx(one,zero,KIND=qp),f(is,1),ldf)

                       call la_wgemm('N','C',mb,js - 1,nb,cmplx(one,zero,KIND=qp),f(is, &
                        js),ldf,e(1,js),lde,cmplx(one,zero,KIND=qp),f(is,1),ldf)

                    end if
                    if (i < p) then
                       call la_wgemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=qp),a( &
                       is,ie + 1),lda,c(is,js),ldc,cmplx(one,zero,KIND=qp),c(ie + 1,js), &
                                 ldc)
                       call la_wgemm('C','N',m - ie,nb,mb,cmplx(-one,zero,KIND=qp),d( &
                       is,ie + 1),ldd,f(is,js),ldf,cmplx(one,zero,KIND=qp),c(ie + 1,js), &
                                 ldc)
                    end if
                 end do loop_200
              end do loop_210
           end if
           work(1) = lwmin
           return
     end subroutine la_wtgsyl
#endif

     !> CTGSEN: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A, B) (in terms of an unitary equivalence trans-
     !> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the pair (A,B). The leading
     !> columns of Q and Z form unitary bases of the corresponding left and
     !> right eigenspaces (deflating subspaces). (A, B) must be in
     !> generalized Schur canonical form, that is, A and B are both upper
     !> triangular.
     !> CTGSEN also computes the generalized eigenvalues
     !> w(j)= ALPHA(j) / BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, the routine computes estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_ctgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alpha,beta,q, &
               ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: dif(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           complex(sp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,ks,liwmin,lwmin,mn2,n1,n2
           real(sp) :: dscale,dsum,rdscal,safmin
           complex(sp) :: temp1,temp2
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,conjg,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -13
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('CTGSEN',-info)
              return
           end if
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
              if (k < n) then
                 if (select(k)) m = m + 1
              else
                 if (select(n)) m = m + 1
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,2*m*(n - m))
              liwmin = max(1,n + 2)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 2)
           else
              lwmin = 1
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -21
           else if (liwork < liwmin .and. .not. lquery) then
              info = -23
           end if
           if (info /= 0) then
              call la_xerbla('CTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_classq(n,a(1,i),1,dscale,dsum)
                    call la_classq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 70
           end if
           ! get machine constant
           safmin = la_slamch('S')
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           do k = 1,n
              swap = select(k)
              if (swap) then
                 ks = ks + 1
                 ! swap the k-th block to position ks. compute unitary q
                 ! and z that will swap adjacent diagonal blocks in (a, b).
                 if (k /= ks) call la_ctgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,k, &
                            ks,ierr)
                 if (ierr > 0) then
                    ! swap is rejected: exit.
                    info = 1
                    if (wantp) then
                       pl = zero
                       pr = zero
                    end if
                    if (wantd) then
                       dif(1) = zero
                       dif(2) = zero
                    end if
                    go to 70
                 end if
              end if
           end do
           if (wantp) then
              ! solve generalized sylvester equation for r and l:
                         ! a11 * r - l * a22 = a12
                         ! b11 * r - l * b22 = b12
              n1 = m
              n2 = n - m
              i = n1 + 1
              call la_clacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_clacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              ijb = 0
              call la_ctgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto
              ! left and right eigenspaces
              rdscal = zero
              dsum = one
              call la_classq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_classq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu estimate.
                 call la_ctgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl estimate.
                 call la_ctgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_clacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_clacn2(mn2,work(mn2 + 1),work,dif(1),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_ctgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_ctgsyl('C',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_clacn2(mn2,work(mn2 + 1),work,dif(2),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_ctgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_ctgsyl('C',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           ! if b(k,k) is complex, make it real and positive (normalization
           ! of the generalized schur form) and store the generalized
           ! eigenvalues of reordered pair (a, b)
           do k = 1,n
              dscale = abs(b(k,k))
              if (dscale > safmin) then
                 temp1 = conjg(b(k,k)/dscale)
                 temp2 = b(k,k)/dscale
                 b(k,k) = dscale
                 call la_cscal(n - k,temp1,b(k,k + 1),ldb)
                 call la_cscal(n - k + 1,temp1,a(k,k),lda)
                 if (wantq) call la_cscal(n,temp2,q(1,k),1)
              else
                 b(k,k) = cmplx(zero,zero,KIND=sp)
              end if
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
           end do
           70 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_ctgsen
     !> ZTGSEN: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A, B) (in terms of an unitary equivalence trans-
     !> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the pair (A,B). The leading
     !> columns of Q and Z form unitary bases of the corresponding left and
     !> right eigenspaces (deflating subspaces). (A, B) must be in
     !> generalized Schur canonical form, that is, A and B are both upper
     !> triangular.
     !> ZTGSEN also computes the generalized eigenvalues
     !> w(j)= ALPHA(j) / BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, the routine computes estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_ztgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alpha,beta,q, &
               ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: dif(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           complex(dp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,ks,liwmin,lwmin,mn2,n1,n2
           real(dp) :: dscale,dsum,rdscal,safmin
           complex(dp) :: temp1,temp2
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,conjg,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -13
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('ZTGSEN',-info)
              return
           end if
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
              if (k < n) then
                 if (select(k)) m = m + 1
              else
                 if (select(n)) m = m + 1
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,2*m*(n - m))
              liwmin = max(1,n + 2)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 2)
           else
              lwmin = 1
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -21
           else if (liwork < liwmin .and. .not. lquery) then
              info = -23
           end if
           if (info /= 0) then
              call la_xerbla('ZTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_zlassq(n,a(1,i),1,dscale,dsum)
                    call la_zlassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 70
           end if
           ! get machine constant
           safmin = la_dlamch('S')
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           do k = 1,n
              swap = select(k)
              if (swap) then
                 ks = ks + 1
                 ! swap the k-th block to position ks. compute unitary q
                 ! and z that will swap adjacent diagonal blocks in (a, b).
                 if (k /= ks) call la_ztgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,k, &
                            ks,ierr)
                 if (ierr > 0) then
                    ! swap is rejected: exit.
                    info = 1
                    if (wantp) then
                       pl = zero
                       pr = zero
                    end if
                    if (wantd) then
                       dif(1) = zero
                       dif(2) = zero
                    end if
                    go to 70
                 end if
              end if
           end do
           if (wantp) then
              ! solve generalized sylvester equation for r and l:
                         ! a11 * r - l * a22 = a12
                         ! b11 * r - l * b22 = b12
              n1 = m
              n2 = n - m
              i = n1 + 1
              call la_zlacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_zlacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              ijb = 0
              call la_ztgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto
              ! left and right eigenspaces
              rdscal = zero
              dsum = one
              call la_zlassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_zlassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu estimate.
                 call la_ztgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl estimate.
                 call la_ztgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_zlacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_zlacn2(mn2,work(mn2 + 1),work,dif(1),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_ztgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_ztgsyl('C',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_zlacn2(mn2,work(mn2 + 1),work,dif(2),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_ztgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_ztgsyl('C',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           ! if b(k,k) is complex, make it real and positive (normalization
           ! of the generalized schur form) and store the generalized
           ! eigenvalues of reordered pair (a, b)
           do k = 1,n
              dscale = abs(b(k,k))
              if (dscale > safmin) then
                 temp1 = conjg(b(k,k)/dscale)
                 temp2 = b(k,k)/dscale
                 b(k,k) = dscale
                 call la_zscal(n - k,temp1,b(k,k + 1),ldb)
                 call la_zscal(n - k + 1,temp1,a(k,k),lda)
                 if (wantq) call la_zscal(n,temp2,q(1,k),1)
              else
                 b(k,k) = cmplx(zero,zero,KIND=dp)
              end if
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
           end do
           70 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_ztgsen
#ifdef LA_WITH_XDP
     !> YTGSEN: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A, B) (in terms of an unitary equivalence trans-
     !> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the pair (A,B). The leading
     !> columns of Q and Z form unitary bases of the corresponding left and
     !> right eigenspaces (deflating subspaces). (A, B) must be in
     !> generalized Schur canonical form, that is, A and B are both upper
     !> triangular.
     !> YTGSEN also computes the generalized eigenvalues
     !> w(j)= ALPHA(j) / BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, the routine computes estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_ytgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alpha,beta,q, &
               ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(xdp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(out) :: dif(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           complex(xdp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,ks,liwmin,lwmin,mn2,n1,n2
           real(xdp) :: dscale,dsum,rdscal,safmin
           complex(xdp) :: temp1,temp2
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,conjg,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -13
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('YTGSEN',-info)
              return
           end if
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
              if (k < n) then
                 if (select(k)) m = m + 1
              else
                 if (select(n)) m = m + 1
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,2*m*(n - m))
              liwmin = max(1,n + 2)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 2)
           else
              lwmin = 1
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -21
           else if (liwork < liwmin .and. .not. lquery) then
              info = -23
           end if
           if (info /= 0) then
              call la_xerbla('YTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_ylassq(n,a(1,i),1,dscale,dsum)
                    call la_ylassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 70
           end if
           ! get machine constant
           safmin = la_xlamch('S')
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           do k = 1,n
              swap = select(k)
              if (swap) then
                 ks = ks + 1
                 ! swap the k-th block to position ks. compute unitary q
                 ! and z that will swap adjacent diagonal blocks in (a, b).
                 if (k /= ks) call la_ytgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,k, &
                            ks,ierr)
                 if (ierr > 0) then
                    ! swap is rejected: exit.
                    info = 1
                    if (wantp) then
                       pl = zero
                       pr = zero
                    end if
                    if (wantd) then
                       dif(1) = zero
                       dif(2) = zero
                    end if
                    go to 70
                 end if
              end if
           end do
           if (wantp) then
              ! solve generalized sylvester equation for r and l:
                         ! a11 * r - l * a22 = a12
                         ! b11 * r - l * b22 = b12
              n1 = m
              n2 = n - m
              i = n1 + 1
              call la_ylacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_ylacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              ijb = 0
              call la_ytgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto
              ! left and right eigenspaces
              rdscal = zero
              dsum = one
              call la_ylassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_ylassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu estimate.
                 call la_ytgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl estimate.
                 call la_ytgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_ylacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_ylacn2(mn2,work(mn2 + 1),work,dif(1),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_ytgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_ytgsyl('C',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_ylacn2(mn2,work(mn2 + 1),work,dif(2),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_ytgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_ytgsyl('C',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           ! if b(k,k) is complex, make it real and positive (normalization
           ! of the generalized schur form) and store the generalized
           ! eigenvalues of reordered pair (a, b)
           do k = 1,n
              dscale = abs(b(k,k))
              if (dscale > safmin) then
                 temp1 = conjg(b(k,k)/dscale)
                 temp2 = b(k,k)/dscale
                 b(k,k) = dscale
                 call la_yscal(n - k,temp1,b(k,k + 1),ldb)
                 call la_yscal(n - k + 1,temp1,a(k,k),lda)
                 if (wantq) call la_yscal(n,temp2,q(1,k),1)
              else
                 b(k,k) = cmplx(zero,zero,KIND=xdp)
              end if
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
           end do
           70 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_ytgsen
#endif
#ifdef LA_WITH_QP
     !> WTGSEN: reorders the generalized Schur decomposition of a complex
     !> matrix pair (A, B) (in terms of an unitary equivalence trans-
     !> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
     !> appears in the leading diagonal blocks of the pair (A,B). The leading
     !> columns of Q and Z form unitary bases of the corresponding left and
     !> right eigenspaces (deflating subspaces). (A, B) must be in
     !> generalized Schur canonical form, that is, A and B are both upper
     !> triangular.
     !> WTGSEN also computes the generalized eigenvalues
     !> w(j)= ALPHA(j) / BETA(j)
     !> of the reordered matrix pair (A, B).
     !> Optionally, the routine computes estimates of reciprocal condition
     !> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
     !> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
     !> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
     !> the selected cluster and the eigenvalues outside the cluster, resp.,
     !> and norms of "projections" onto left and right eigenspaces w.r.t.
     !> the selected cluster in the (1,1)-block.

     pure subroutine la_wtgsen(ijob,wantq,wantz,select,n,a,lda,b,ldb,alpha,beta,q, &
               ldq,z,ldz,m,pl,pr,dif,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq,wantz
           integer(ilp),intent(in) :: ijob,lda,ldb,ldq,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(out) :: pl,pr
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: dif(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),z(ldz,*)
           complex(qp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,swap,wantd,wantd1,wantd2,wantp
           integer(ilp) :: i,ierr,ijb,k,kase,ks,liwmin,lwmin,mn2,n1,n2
           real(qp) :: dscale,dsum,rdscal,safmin
           complex(qp) :: temp1,temp2
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,conjg,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
           if (ijob < 0 .or. ijob > 5) then
              info = -1
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -13
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('WTGSEN',-info)
              return
           end if
           ierr = 0
           wantp = ijob == 1 .or. ijob >= 4
           wantd1 = ijob == 2 .or. ijob == 4
           wantd2 = ijob == 3 .or. ijob == 5
           wantd = wantd1 .or. wantd2
           ! set m to the dimension of the specified pair of deflating
           ! subspaces.
           m = 0
           if (.not. lquery .or. ijob /= 0) then
           do k = 1,n
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
              if (k < n) then
                 if (select(k)) m = m + 1
              else
                 if (select(n)) m = m + 1
              end if
           end do
           end if
           if (ijob == 1 .or. ijob == 2 .or. ijob == 4) then
              lwmin = max(1,2*m*(n - m))
              liwmin = max(1,n + 2)
           else if (ijob == 3 .or. ijob == 5) then
              lwmin = max(1,4*m*(n - m))
              liwmin = max(1,2*m*(n - m),n + 2)
           else
              lwmin = 1
              liwmin = 1
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           if (lwork < lwmin .and. .not. lquery) then
              info = -21
           else if (liwork < liwmin .and. .not. lquery) then
              info = -23
           end if
           if (info /= 0) then
              call la_xerbla('WTGSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wantp) then
                 pl = one
                 pr = one
              end if
              if (wantd) then
                 dscale = zero
                 dsum = one
                 do i = 1,n
                    call la_wlassq(n,a(1,i),1,dscale,dsum)
                    call la_wlassq(n,b(1,i),1,dscale,dsum)
                 end do
                 dif(1) = dscale*sqrt(dsum)
                 dif(2) = dif(1)
              end if
              go to 70
           end if
           ! get machine constant
           safmin = la_qlamch('S')
           ! collect the selected blocks at the top-left corner of (a, b).
           ks = 0
           do k = 1,n
              swap = select(k)
              if (swap) then
                 ks = ks + 1
                 ! swap the k-th block to position ks. compute unitary q
                 ! and z that will swap adjacent diagonal blocks in (a, b).
                 if (k /= ks) call la_wtgexc(wantq,wantz,n,a,lda,b,ldb,q,ldq,z,ldz,k, &
                            ks,ierr)
                 if (ierr > 0) then
                    ! swap is rejected: exit.
                    info = 1
                    if (wantp) then
                       pl = zero
                       pr = zero
                    end if
                    if (wantd) then
                       dif(1) = zero
                       dif(2) = zero
                    end if
                    go to 70
                 end if
              end if
           end do
           if (wantp) then
              ! solve generalized sylvester equation for r and l:
                         ! a11 * r - l * a22 = a12
                         ! b11 * r - l * b22 = b12
              n1 = m
              n2 = n - m
              i = n1 + 1
              call la_wlacpy('FULL',n1,n2,a(1,i),lda,work,n1)
              call la_wlacpy('FULL',n1,n2,b(1,i),ldb,work(n1*n2 + 1),n1)
              ijb = 0
              call la_wtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b(i, &
               i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - 2*n1*n2, &
                         iwork,ierr)
              ! estimate the reciprocal of norms of "projections" onto
              ! left and right eigenspaces
              rdscal = zero
              dsum = one
              call la_wlassq(n1*n2,work,1,rdscal,dsum)
              pl = rdscal*sqrt(dsum)
              if (pl == zero) then
                 pl = one
              else
                 pl = dscale/(sqrt(dscale*dscale/pl + pl)*sqrt(pl))
              end if
              rdscal = zero
              dsum = one
              call la_wlassq(n1*n2,work(n1*n2 + 1),1,rdscal,dsum)
              pr = rdscal*sqrt(dsum)
              if (pr == zero) then
                 pr = one
              else
                 pr = dscale/(sqrt(dscale*dscale/pr + pr)*sqrt(pr))
              end if
           end if
           if (wantd) then
              ! compute estimates difu and difl.
              if (wantd1) then
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = idifjb
                 ! frobenius norm-based difu estimate.
                 call la_wtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b,ldb,b( &
                  i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
                 ! frobenius norm-based difl estimate.
                 call la_wtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b(i,i), &
                  ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1),lwork - &
                            2*n1*n2,iwork,ierr)
              else
                 ! compute 1-norm-based estimates of difu and difl using
                 ! reversed communication with la_wlacn2. in each step a
                 ! generalized sylvester equation or a transposed variant
                 ! is solved.
                 kase = 0
                 n1 = m
                 n2 = n - m
                 i = n1 + 1
                 ijb = 0
                 mn2 = 2*n1*n2
                 ! 1-norm-based estimate of difu.
                 40 continue
                 call la_wlacn2(mn2,work(mn2 + 1),work,dif(1),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_wtgsyl('N',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_wtgsyl('C',ijb,n1,n2,a,lda,a(i,i),lda,work,n1,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n1,dscale,dif(1),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 40
                 end if
                 dif(1) = dscale/dif(1)
                 ! 1-norm-based estimate of difl.
                 50 continue
                 call la_wlacn2(mn2,work(mn2 + 1),work,dif(2),kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve generalized sylvester equation
                       call la_wtgsyl('N',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b( &
                       i,i),ldb,b,ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    else
                       ! solve the transposed variant.
                       call la_wtgsyl('C',ijb,n2,n1,a(i,i),lda,a,lda,work,n2,b, &
                       ldb,b(i,i),ldb,work(n1*n2 + 1),n2,dscale,dif(2),work(n1*n2*2 + 1) &
                                 ,lwork - 2*n1*n2,iwork,ierr)
                    end if
                    go to 50
                 end if
                 dif(2) = dscale/dif(2)
              end if
           end if
           ! if b(k,k) is complex, make it real and positive (normalization
           ! of the generalized schur form) and store the generalized
           ! eigenvalues of reordered pair (a, b)
           do k = 1,n
              dscale = abs(b(k,k))
              if (dscale > safmin) then
                 temp1 = conjg(b(k,k)/dscale)
                 temp2 = b(k,k)/dscale
                 b(k,k) = dscale
                 call la_wscal(n - k,temp1,b(k,k + 1),ldb)
                 call la_wscal(n - k + 1,temp1,a(k,k),lda)
                 if (wantq) call la_wscal(n,temp2,q(1,k),1)
              else
                 b(k,k) = cmplx(zero,zero,KIND=qp)
              end if
              alpha(k) = a(k,k)
              beta(k) = b(k,k)
           end do
           70 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_wtgsen
#endif

     !> CTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B).
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.

     pure subroutine la_ctgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: dif(*),s(*)
           complex(sp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,k,ks,lwmin,n1,n2
           real(sp) :: bignum,cond,eps,lnrm,rnrm,scale,smlnum
           complex(sp) :: yhax,yhbx
           ! Local Arrays
           complex(sp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,max
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 do k = 1,n
                    if (select(k)) m = m + 1
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*n
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ks = 0
           loop_20: do k = 1,n
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (.not. select(k)) cycle loop_20
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 rnrm = la_scnrm2(n,vr(1,ks),1)
                 lnrm = la_scnrm2(n,vl(1,ks),1)
                 call la_cgemv('N',n,n,cmplx(one,zero,KIND=sp),a,lda,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=sp),work,1)
                 yhax = la_cdotc(n,work,1,vl(1,ks),1)
                 call la_cgemv('N',n,n,cmplx(one,zero,KIND=sp),b,ldb,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=sp),work,1)
                 yhbx = la_cdotc(n,work,1,vl(1,ks),1)
                 cond = la_slapy2(abs(yhax),abs(yhbx))
                 if (cond == zero) then
                    s(ks) = -one
                 else
                    s(ks) = cond/(rnrm*lnrm)
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_slapy2(abs(a(1,1)),abs(b(1,1)))
                 else
                    ! estimate the reciprocal condition number of the k-th
                    ! eigenvectors.
                    ! copy the matrix (a, b) to the array work and move the
                    ! (k,k)th pair to the (1,1) position.
                    call la_clacpy('FULL',n,n,a,lda,work,n)
                    call la_clacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                    ifst = k
                    ilst = 1
                    call la_ctgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                              dummy1,1,ifst,ilst,ierr)
                    if (ierr > 0) then
                       ! ill-conditioned problem - swap rejected.
                       dif(ks) = zero
                    else
                       ! reordering successful, solve generalized sylvester
                       ! equation for r and l,
                                  ! a22 * r - l * a11 = a12
                                  ! b22 * r - l * b11 = b12,
                       ! and compute estimate of difl[(a11,b11), (a22, b22)].
                       n1 = 1
                       n2 = n - n1
                       i = n*n + 1
                       call la_ctgsyl('N',idifjb,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),dummy,1,iwork,ierr)
                    end if
                 end if
              end if
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_ctgsna
     !> ZTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B).
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.

     pure subroutine la_ztgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: dif(*),s(*)
           complex(dp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,k,ks,lwmin,n1,n2
           real(dp) :: bignum,cond,eps,lnrm,rnrm,scale,smlnum
           complex(dp) :: yhax,yhbx
           ! Local Arrays
           complex(dp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,max
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 do k = 1,n
                    if (select(k)) m = m + 1
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*n
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ks = 0
           loop_20: do k = 1,n
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (.not. select(k)) cycle loop_20
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 rnrm = la_dznrm2(n,vr(1,ks),1)
                 lnrm = la_dznrm2(n,vl(1,ks),1)
                 call la_zgemv('N',n,n,cmplx(one,zero,KIND=dp),a,lda,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=dp),work,1)
                 yhax = la_zdotc(n,work,1,vl(1,ks),1)
                 call la_zgemv('N',n,n,cmplx(one,zero,KIND=dp),b,ldb,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=dp),work,1)
                 yhbx = la_zdotc(n,work,1,vl(1,ks),1)
                 cond = la_dlapy2(abs(yhax),abs(yhbx))
                 if (cond == zero) then
                    s(ks) = -one
                 else
                    s(ks) = cond/(rnrm*lnrm)
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_dlapy2(abs(a(1,1)),abs(b(1,1)))
                 else
                    ! estimate the reciprocal condition number of the k-th
                    ! eigenvectors.
                    ! copy the matrix (a, b) to the array work and move the
                    ! (k,k)th pair to the (1,1) position.
                    call la_zlacpy('FULL',n,n,a,lda,work,n)
                    call la_zlacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                    ifst = k
                    ilst = 1
                    call la_ztgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                              dummy1,1,ifst,ilst,ierr)
                    if (ierr > 0) then
                       ! ill-conditioned problem - swap rejected.
                       dif(ks) = zero
                    else
                       ! reordering successful, solve generalized sylvester
                       ! equation for r and l,
                                  ! a22 * r - l * a11 = a12
                                  ! b22 * r - l * b11 = b12,
                       ! and compute estimate of difl[(a11,b11), (a22, b22)].
                       n1 = 1
                       n2 = n - n1
                       i = n*n + 1
                       call la_ztgsyl('N',idifjb,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),dummy,1,iwork,ierr)
                    end if
                 end if
              end if
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_ztgsna
#ifdef LA_WITH_XDP
     !> YTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B).
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.

     pure subroutine la_ytgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(out) :: dif(*),s(*)
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,k,ks,lwmin,n1,n2
           real(xdp) :: bignum,cond,eps,lnrm,rnrm,scale,smlnum
           complex(xdp) :: yhax,yhbx
           ! Local Arrays
           complex(xdp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,max
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 do k = 1,n
                    if (select(k)) m = m + 1
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*n
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ks = 0
           loop_20: do k = 1,n
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (.not. select(k)) cycle loop_20
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 rnrm = la_xynrm2(n,vr(1,ks),1)
                 lnrm = la_xynrm2(n,vl(1,ks),1)
                 call la_ygemv('N',n,n,cmplx(one,zero,KIND=xdp),a,lda,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=xdp),work,1)
                 yhax = la_ydotc(n,work,1,vl(1,ks),1)
                 call la_ygemv('N',n,n,cmplx(one,zero,KIND=xdp),b,ldb,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=xdp),work,1)
                 yhbx = la_ydotc(n,work,1,vl(1,ks),1)
                 cond = la_xlapy2(abs(yhax),abs(yhbx))
                 if (cond == zero) then
                    s(ks) = -one
                 else
                    s(ks) = cond/(rnrm*lnrm)
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_xlapy2(abs(a(1,1)),abs(b(1,1)))
                 else
                    ! estimate the reciprocal condition number of the k-th
                    ! eigenvectors.
                    ! copy the matrix (a, b) to the array work and move the
                    ! (k,k)th pair to the (1,1) position.
                    call la_ylacpy('FULL',n,n,a,lda,work,n)
                    call la_ylacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                    ifst = k
                    ilst = 1
                    call la_ytgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                              dummy1,1,ifst,ilst,ierr)
                    if (ierr > 0) then
                       ! ill-conditioned problem - swap rejected.
                       dif(ks) = zero
                    else
                       ! reordering successful, solve generalized sylvester
                       ! equation for r and l,
                                  ! a22 * r - l * a11 = a12
                                  ! b22 * r - l * b11 = b12,
                       ! and compute estimate of difl[(a11,b11), (a22, b22)].
                       n1 = 1
                       n2 = n - n1
                       i = n*n + 1
                       call la_ytgsyl('N',idifjb,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),dummy,1,iwork,ierr)
                    end if
                 end if
              end if
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_ytgsna
#endif
#ifdef LA_WITH_QP
     !> WTGSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or eigenvectors of a matrix pair (A, B).
     !> (A, B) must be in generalized Schur canonical form, that is, A and
     !> B are both upper triangular.

     pure subroutine la_wtgsna(job,howmny,select,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,s, &
               dif,mm,m,work,lwork,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: dif(*),s(*)
           complex(qp),intent(in) :: a(lda,*),b(ldb,*),vl(ldvl,*),vr(ldvr,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: idifjb = 3

           ! Local Scalars
           logical(lk) :: lquery,somcon,wantbh,wantdf,wants
           integer(ilp) :: i,ierr,ifst,ilst,k,ks,lwmin,n1,n2
           real(qp) :: bignum,cond,eps,lnrm,rnrm,scale,smlnum
           complex(qp) :: yhax,yhbx
           ! Local Arrays
           complex(qp) :: dummy(1),dummy1(1)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,max
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantdf = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           lquery = (lwork == -1)
           if (.not. wants .and. .not. wantdf) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           else if (wants .and. ldvl < n) then
              info = -10
           else if (wants .and. ldvr < n) then
              info = -12
           else
              ! set m to the number of eigenpairs for which condition numbers
              ! are required, and test mm.
              if (somcon) then
                 m = 0
                 do k = 1,n
                    if (select(k)) m = m + 1
                 end do
              else
                 m = n
              end if
              if (n == 0) then
                 lwmin = 1
              else if (la_lsame(job,'V') .or. la_lsame(job,'B')) then
                 lwmin = 2*n*n
              else
                 lwmin = n
              end if
              work(1) = lwmin
              if (mm < m) then
                 info = -15
              else if (lwork < lwmin .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WTGSNA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ks = 0
           loop_20: do k = 1,n
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (.not. select(k)) cycle loop_20
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 rnrm = la_qwnrm2(n,vr(1,ks),1)
                 lnrm = la_qwnrm2(n,vl(1,ks),1)
                 call la_wgemv('N',n,n,cmplx(one,zero,KIND=qp),a,lda,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=qp),work,1)
                 yhax = la_wdotc(n,work,1,vl(1,ks),1)
                 call la_wgemv('N',n,n,cmplx(one,zero,KIND=qp),b,ldb,vr(1,ks),1, &
                           cmplx(zero,zero,KIND=qp),work,1)
                 yhbx = la_wdotc(n,work,1,vl(1,ks),1)
                 cond = la_qlapy2(abs(yhax),abs(yhbx))
                 if (cond == zero) then
                    s(ks) = -one
                 else
                    s(ks) = cond/(rnrm*lnrm)
                 end if
              end if
              if (wantdf) then
                 if (n == 1) then
                    dif(ks) = la_qlapy2(abs(a(1,1)),abs(b(1,1)))
                 else
                    ! estimate the reciprocal condition number of the k-th
                    ! eigenvectors.
                    ! copy the matrix (a, b) to the array work and move the
                    ! (k,k)th pair to the (1,1) position.
                    call la_wlacpy('FULL',n,n,a,lda,work,n)
                    call la_wlacpy('FULL',n,n,b,ldb,work(n*n + 1),n)
                    ifst = k
                    ilst = 1
                    call la_wtgexc(.false.,.false.,n,work,n,work(n*n + 1),n,dummy,1, &
                              dummy1,1,ifst,ilst,ierr)
                    if (ierr > 0) then
                       ! ill-conditioned problem - swap rejected.
                       dif(ks) = zero
                    else
                       ! reordering successful, solve generalized sylvester
                       ! equation for r and l,
                                  ! a22 * r - l * a11 = a12
                                  ! b22 * r - l * b11 = b12,
                       ! and compute estimate of difl[(a11,b11), (a22, b22)].
                       n1 = 1
                       n2 = n - n1
                       i = n*n + 1
                       call la_wtgsyl('N',idifjb,n2,n1,work(n*n1 + n1 + 1),n,work,n, &
                       work(n1 + 1),n,work(n*n1 + n1 + i),n,work(i),n,work(n1 + i),n,scale, &
                                 dif(ks),dummy,1,iwork,ierr)
                    end if
                 end if
              end if
           end do loop_20
           work(1) = lwmin
           return
     end subroutine la_wtgsna
#endif

end module la_lapack_eigv_comp2
