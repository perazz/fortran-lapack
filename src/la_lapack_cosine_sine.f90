!> Cosine-sine decomposition: bidiagonal block form, simultaneous bidiagonalization, row and column permutations
module la_lapack_cosine_sine
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_givens_jacobi_rot
     use la_lapack_householder_reflectors
     use la_lapack_orthogonal_factors_ql
     use la_lapack_orthogonal_factors_qr
     use la_lapack_svd_comp2
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slapmr
     public :: la_slapmt
     public :: la_sorbdb6
     public :: la_sbbcsd
     public :: la_sorbdb
     public :: la_sorbdb5
     public :: la_sorcsd
     public :: la_sorbdb1
     public :: la_sorbdb2
     public :: la_sorbdb3
     public :: la_sorbdb4
     public :: la_sorcsd2by1
     public :: la_dlapmr
     public :: la_dlapmt
     public :: la_dorbdb6
     public :: la_dbbcsd
     public :: la_dorbdb
     public :: la_dorbdb5
     public :: la_dorcsd
     public :: la_dorbdb1
     public :: la_dorbdb2
     public :: la_dorbdb3
     public :: la_dorbdb4
     public :: la_dorcsd2by1
     public :: la_qlapmr
     public :: la_qlapmt
     public :: la_qorbdb6
     public :: la_qbbcsd
     public :: la_qorbdb
     public :: la_qorbdb5
     public :: la_qorcsd
     public :: la_qorbdb1
     public :: la_qorbdb2
     public :: la_qorbdb3
     public :: la_qorbdb4
     public :: la_qorcsd2by1
     public :: la_clapmr
     public :: la_clapmt
     public :: la_cunbdb
     public :: la_cunbdb6
     public :: la_cbbcsd
     public :: la_cunbdb5
     public :: la_cuncsd
     public :: la_cunbdb1
     public :: la_cunbdb2
     public :: la_cunbdb3
     public :: la_cunbdb4
     public :: la_cuncsd2by1
     public :: la_zlapmr
     public :: la_zlapmt
     public :: la_zunbdb
     public :: la_zunbdb6
     public :: la_zbbcsd
     public :: la_zunbdb5
     public :: la_zuncsd
     public :: la_zunbdb1
     public :: la_zunbdb2
     public :: la_zunbdb3
     public :: la_zunbdb4
     public :: la_zuncsd2by1
     public :: la_wlapmr
     public :: la_wlapmt
     public :: la_wunbdb
     public :: la_wunbdb6
     public :: la_wbbcsd
     public :: la_wunbdb5
     public :: la_wuncsd
     public :: la_wunbdb1
     public :: la_wunbdb2
     public :: la_wunbdb3
     public :: la_wunbdb4
     public :: la_wuncsd2by1

     contains

     !> SLAPMR: rearranges the rows of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.

     pure subroutine la_slapmr(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           real(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,in,j,jj
           real(sp) :: temp
           ! Executable Statements
           if (m <= 1) return
           do i = 1,m
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,m
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do jj = 1,n
                    temp = x(j,jj)
                    x(j,jj) = x(in,jj)
                    x(in,jj) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,m
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do jj = 1,n
                    temp = x(i,jj)
                    x(i,jj) = x(j,jj)
                    x(j,jj) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_slapmr
     !> DLAPMR: rearranges the rows of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.

     pure subroutine la_dlapmr(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           real(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,in,j,jj
           real(dp) :: temp
           ! Executable Statements
           if (m <= 1) return
           do i = 1,m
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,m
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do jj = 1,n
                    temp = x(j,jj)
                    x(j,jj) = x(in,jj)
                    x(in,jj) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,m
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do jj = 1,n
                    temp = x(i,jj)
                    x(i,jj) = x(j,jj)
                    x(j,jj) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_dlapmr
     !> QLAPMR: rearranges the rows of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.

     pure subroutine la_qlapmr(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           real(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,in,j,jj
           real(qp) :: temp
           ! Executable Statements
           if (m <= 1) return
           do i = 1,m
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,m
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do jj = 1,n
                    temp = x(j,jj)
                    x(j,jj) = x(in,jj)
                    x(in,jj) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,m
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do jj = 1,n
                    temp = x(i,jj)
                    x(i,jj) = x(j,jj)
                    x(j,jj) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_qlapmr

     !> SLAPMT: rearranges the columns of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.

     pure subroutine la_slapmt(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           real(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ii,j,in
           real(sp) :: temp
           ! Executable Statements
           if (n <= 1) return
           do i = 1,n
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,n
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do ii = 1,m
                    temp = x(ii,j)
                    x(ii,j) = x(ii,in)
                    x(ii,in) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,n
                 if (k(i) > 0) go to 100
                 k(i) = -k(i)
                 j = k(i)
                 80 continue
                 if (j == i) go to 100
                 do ii = 1,m
                    temp = x(ii,i)
                    x(ii,i) = x(ii,j)
                    x(ii,j) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 80
                 100 continue
              end do
           end if
           return
     end subroutine la_slapmt
     !> DLAPMT: rearranges the columns of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.

     pure subroutine la_dlapmt(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           real(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ii,in,j
           real(dp) :: temp
           ! Executable Statements
           if (n <= 1) return
           do i = 1,n
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,n
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do ii = 1,m
                    temp = x(ii,j)
                    x(ii,j) = x(ii,in)
                    x(ii,in) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,n
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do ii = 1,m
                    temp = x(ii,i)
                    x(ii,i) = x(ii,j)
                    x(ii,j) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_dlapmt
     !> QLAPMT: rearranges the columns of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.

     pure subroutine la_qlapmt(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           real(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ii,in,j
           real(qp) :: temp
           ! Executable Statements
           if (n <= 1) return
           do i = 1,n
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,n
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do ii = 1,m
                    temp = x(ii,j)
                    x(ii,j) = x(ii,in)
                    x(ii,in) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,n
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do ii = 1,m
                    temp = x(ii,i)
                    x(ii,i) = x(ii,j)
                    x(ii,j) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_qlapmt

     !> SORBDB6: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then the zero vector is returned.

     pure subroutine la_sorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           real(sp),intent(out) :: work(*)
           real(sp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: alphasq = 0.01_sp
           real(sp),parameter :: realone = 1.0_sp
           real(sp),parameter :: realzero = 0.0_sp

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: normsq1,normsq2,scl1,scl2,ssq1,ssq2
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('SORBDB6',-info)
              return
           end if
           ! first, project x onto the orthogonal complement of q's column
           ! space
           scl1 = realzero
           ssq1 = realone
           call la_slassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_slassq(m2,x2,incx2,scl2,ssq2)
           normsq1 = scl1**2*ssq1 + scl2**2*ssq2
           if (m1 == 0) then
              do i = 1,n
                 work(i) = zero
              end do
           else
              call la_sgemv('C',m1,n,one,q1,ldq1,x1,incx1,zero,work,1)
           end if
           call la_sgemv('C',m2,n,one,q2,ldq2,x2,incx2,one,work,1)
           call la_sgemv('N',m1,n,negone,q1,ldq1,work,1,one,x1,incx1)
           call la_sgemv('N',m2,n,negone,q2,ldq2,work,1,one,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_slassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_slassq(m2,x2,incx2,scl2,ssq2)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if projection is sufficiently large in norm, then stop.
           ! if projection is zero, then stop.
           ! otherwise, project again.
           if (normsq2 >= alphasq*normsq1) then
              return
           end if
           if (normsq2 == zero) then
              return
           end if
           normsq1 = normsq2
           do i = 1,n
              work(i) = zero
           end do
           if (m1 == 0) then
              do i = 1,n
                 work(i) = zero
              end do
           else
              call la_sgemv('C',m1,n,one,q1,ldq1,x1,incx1,zero,work,1)
           end if
           call la_sgemv('C',m2,n,one,q2,ldq2,x2,incx2,one,work,1)
           call la_sgemv('N',m1,n,negone,q1,ldq1,work,1,one,x1,incx1)
           call la_sgemv('N',m2,n,negone,q2,ldq2,work,1,one,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_slassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_slassq(m1,x1,incx1,scl1,ssq1)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if second projection is sufficiently large in norm, then do
           ! nothing more. alternatively, if it shrunk significantly, then
           ! truncate it to zero.
           if (normsq2 < alphasq*normsq1) then
              do i = 1,m1
                 x1(i) = zero
              end do
              do i = 1,m2
                 x2(i) = zero
              end do
           end if
           return
     end subroutine la_sorbdb6
     !> DORBDB6: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then the zero vector is returned.

     pure subroutine la_dorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           real(dp),intent(out) :: work(*)
           real(dp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: alphasq = 0.01_dp
           real(dp),parameter :: realone = 1.0_dp
           real(dp),parameter :: realzero = 0.0_dp

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: normsq1,normsq2,scl1,scl2,ssq1,ssq2
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('DORBDB6',-info)
              return
           end if
           ! first, project x onto the orthogonal complement of q's column
           ! space
           scl1 = realzero
           ssq1 = realone
           call la_dlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_dlassq(m2,x2,incx2,scl2,ssq2)
           normsq1 = scl1**2*ssq1 + scl2**2*ssq2
           if (m1 == 0) then
              do i = 1,n
                 work(i) = zero
              end do
           else
              call la_dgemv('C',m1,n,one,q1,ldq1,x1,incx1,zero,work,1)
           end if
           call la_dgemv('C',m2,n,one,q2,ldq2,x2,incx2,one,work,1)
           call la_dgemv('N',m1,n,negone,q1,ldq1,work,1,one,x1,incx1)
           call la_dgemv('N',m2,n,negone,q2,ldq2,work,1,one,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_dlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_dlassq(m2,x2,incx2,scl2,ssq2)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if projection is sufficiently large in norm, then stop.
           ! if projection is zero, then stop.
           ! otherwise, project again.
           if (normsq2 >= alphasq*normsq1) then
              return
           end if
           if (normsq2 == zero) then
              return
           end if
           normsq1 = normsq2
           do i = 1,n
              work(i) = zero
           end do
           if (m1 == 0) then
              do i = 1,n
                 work(i) = zero
              end do
           else
              call la_dgemv('C',m1,n,one,q1,ldq1,x1,incx1,zero,work,1)
           end if
           call la_dgemv('C',m2,n,one,q2,ldq2,x2,incx2,one,work,1)
           call la_dgemv('N',m1,n,negone,q1,ldq1,work,1,one,x1,incx1)
           call la_dgemv('N',m2,n,negone,q2,ldq2,work,1,one,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_dlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_dlassq(m1,x1,incx1,scl1,ssq1)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if second projection is sufficiently large in norm, then do
           ! nothing more. alternatively, if it shrunk significantly, then
           ! truncate it to zero.
           if (normsq2 < alphasq*normsq1) then
              do i = 1,m1
                 x1(i) = zero
              end do
              do i = 1,m2
                 x2(i) = zero
              end do
           end if
           return
     end subroutine la_dorbdb6
     !> QORBDB6: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then the zero vector is returned.

     pure subroutine la_qorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           real(qp),intent(out) :: work(*)
           real(qp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: alphasq = 0.01_qp
           real(qp),parameter :: realone = 1.0_qp
           real(qp),parameter :: realzero = 0.0_qp

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: normsq1,normsq2,scl1,scl2,ssq1,ssq2
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('QORBDB6',-info)
              return
           end if
           ! first, project x onto the orthogonal complement of q's column
           ! space
           scl1 = realzero
           ssq1 = realone
           call la_qlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_qlassq(m2,x2,incx2,scl2,ssq2)
           normsq1 = scl1**2*ssq1 + scl2**2*ssq2
           if (m1 == 0) then
              do i = 1,n
                 work(i) = zero
              end do
           else
              call la_qgemv('C',m1,n,one,q1,ldq1,x1,incx1,zero,work,1)
           end if
           call la_qgemv('C',m2,n,one,q2,ldq2,x2,incx2,one,work,1)
           call la_qgemv('N',m1,n,negone,q1,ldq1,work,1,one,x1,incx1)
           call la_qgemv('N',m2,n,negone,q2,ldq2,work,1,one,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_qlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_qlassq(m2,x2,incx2,scl2,ssq2)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if projection is sufficiently large in norm, then stop.
           ! if projection is zero, then stop.
           ! otherwise, project again.
           if (normsq2 >= alphasq*normsq1) then
              return
           end if
           if (normsq2 == zero) then
              return
           end if
           normsq1 = normsq2
           do i = 1,n
              work(i) = zero
           end do
           if (m1 == 0) then
              do i = 1,n
                 work(i) = zero
              end do
           else
              call la_qgemv('C',m1,n,one,q1,ldq1,x1,incx1,zero,work,1)
           end if
           call la_qgemv('C',m2,n,one,q2,ldq2,x2,incx2,one,work,1)
           call la_qgemv('N',m1,n,negone,q1,ldq1,work,1,one,x1,incx1)
           call la_qgemv('N',m2,n,negone,q2,ldq2,work,1,one,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_qlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_qlassq(m1,x1,incx1,scl1,ssq1)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if second projection is sufficiently large in norm, then do
           ! nothing more. alternatively, if it shrunk significantly, then
           ! truncate it to zero.
           if (normsq2 < alphasq*normsq1) then
              do i = 1,m1
                 x1(i) = zero
              end do
              do i = 1,m2
                 x2(i) = zero
              end do
           end if
           return
     end subroutine la_qorbdb6

     !> SBBCSD: computes the CS decomposition of an orthogonal matrix in
     !> bidiagonal-block form,
     !> [ B11 | B12 0  0 ]
     !> [  0  |  0 -I  0 ]
     !> X = [----------------]
     !> [ B21 | B22 0  0 ]
     !> [  0  |  0  0  I ]
     !> [  C | -S  0  0 ]
     !> [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T
     !> = [---------] [---------------] [---------]   .
     !> [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
     !> [  0 |  0  0  I ]
     !> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
     !> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
     !> transposed and/or permuted. This can be done in constant time using
     !> the TRANS and SIGNS options. See SORCSD for details.)
     !> The bidiagonal matrices B11, B12, B21, and B22 are represented
     !> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
     !> The orthogonal matrices U1, U2, V1T, and V2T are input/output.
     !> The input matrices are pre- or post-multiplied by the appropriate
     !> singular vector matrices.

     pure subroutine la_sbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,phi,u1, &
     ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,b11d,b11e,b12d,b12e,b21d,b21e,b22d,b22e,work, &
               lwork,info)
        use la_constants_sp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,lwork,m,p,q
           ! Array Arguments
           real(sp),intent(out) :: b11d(*),b11e(*),b12d(*),b12e(*),b21d(*),b21e(*),b22d(*), &
                      b22e(*),work(*)
           real(sp),intent(inout) :: phi(*),theta(*)
           real(sp),intent(inout) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*)
        ! ===================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 6
           real(sp),parameter :: hundred = 100.0_sp
           real(sp),parameter :: meighth = -0.125_sp
           real(sp),parameter :: piover2 = 1.57079632679489661923132169163975144210_sp

           ! Local Scalars
           logical(lk) :: colmajor,lquery,restart11,restart12,restart21,restart22,wantu1, &
                     wantu2,wantv1t,wantv2t
           integer(ilp) :: i,imin,imax,iter,iu1cs,iu1sn,iu2cs,iu2sn,iv1tcs,iv1tsn, &
                     iv2tcs,iv2tsn,j,lworkmin,lworkopt,maxit,mini
           real(sp) :: b11bulge,b12bulge,b21bulge,b22bulge,dummy,eps,mu,nu,r,sigma11, &
                     sigma21,temp,thetamax,thetamin,thresh,tol,tolmul,unfl,x1,x2,y1,y2
           ! Intrinsic Functions
           intrinsic :: abs,atan2,cos,max,min,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           if (m < 0) then
              info = -6
           else if (p < 0 .or. p > m) then
              info = -7
           else if (q < 0 .or. q > m) then
              info = -8
           else if (q > p .or. q > m - p .or. q > m - q) then
              info = -8
           else if (wantu1 .and. ldu1 < p) then
              info = -12
           else if (wantu2 .and. ldu2 < m - p) then
              info = -14
           else if (wantv1t .and. ldv1t < q) then
              info = -16
           else if (wantv2t .and. ldv2t < m - q) then
              info = -18
           end if
           ! quick return if q = 0
           if (info == 0 .and. q == 0) then
              lworkmin = 1
              work(1) = lworkmin
              return
           end if
           ! compute workspace
           if (info == 0) then
              iu1cs = 1
              iu1sn = iu1cs + q
              iu2cs = iu1sn + q
              iu2sn = iu2cs + q
              iv1tcs = iu2sn + q
              iv1tsn = iv1tcs + q
              iv2tcs = iv1tsn + q
              iv2tsn = iv2tcs + q
              lworkopt = iv2tsn + q - 1
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -28
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SBBCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! get machine constants
           eps = la_slamch('EPSILON')
           unfl = la_slamch('SAFE MINIMUM')
           tolmul = max(ten,min(hundred,eps**meighth))
           tol = tolmul*eps
           thresh = max(tol,maxitr*q*q*unfl)
           ! test for negligible sines or cosines
           do i = 1,q
              if (theta(i) < thresh) then
                 theta(i) = zero
              else if (theta(i) > piover2 - thresh) then
                 theta(i) = piover2
              end if
           end do
           do i = 1,q - 1
              if (phi(i) < thresh) then
                 phi(i) = zero
              else if (phi(i) > piover2 - thresh) then
                 phi(i) = piover2
              end if
           end do
           ! initial deflation
           imax = q
           do while (imax > 1)
              if (phi(imax - 1) /= zero) then
                 exit
              end if
              imax = imax - 1
           end do
           imin = imax - 1
           if (imin > 1) then
              do while (phi(imin - 1) /= zero)
                 imin = imin - 1
                 if (imin <= 1) exit
              end do
           end if
           ! initialize iteration counter
           maxit = maxitr*q*q
           iter = 0
           ! begin main iteration loop
           do while (imax > 1)
              ! compute the matrix entries
              b11d(imin) = cos(theta(imin))
              b21d(imin) = -sin(theta(imin))
              do i = imin,imax - 1
                 b11e(i) = -sin(theta(i))*sin(phi(i))
                 b11d(i + 1) = cos(theta(i + 1))*cos(phi(i))
                 b12d(i) = sin(theta(i))*cos(phi(i))
                 b12e(i) = cos(theta(i + 1))*sin(phi(i))
                 b21e(i) = -cos(theta(i))*sin(phi(i))
                 b21d(i + 1) = -sin(theta(i + 1))*cos(phi(i))
                 b22d(i) = cos(theta(i))*cos(phi(i))
                 b22e(i) = -sin(theta(i + 1))*sin(phi(i))
              end do
              b12d(imax) = sin(theta(imax))
              b22d(imax) = cos(theta(imax))
              ! abort if not converging; otherwise, increment iter
              if (iter > maxit) then
                 info = 0
                 do i = 1,q
                    if (phi(i) /= zero) info = info + 1
                 end do
                 return
              end if
              iter = iter + imax - imin
              ! compute shifts
              thetamax = theta(imin)
              thetamin = theta(imin)
              do i = imin + 1,imax
                 if (theta(i) > thetamax) thetamax = theta(i)
                 if (theta(i) < thetamin) thetamin = theta(i)
              end do
              if (thetamax > piover2 - thresh) then
                 ! zero on diagonals of b11 and b22; induce deflation with a
                 ! zero shift
                 mu = zero
                 nu = one
              else if (thetamin < thresh) then
                 ! zero on diagonals of b12 and b22; induce deflation with a
                 ! zero shift
                 mu = one
                 nu = zero
              else
                 ! compute shifts for b11 and b21 and use the lesser
                 call la_slas2(b11d(imax - 1),b11e(imax - 1),b11d(imax),sigma11,dummy)

                 call la_slas2(b21d(imax - 1),b21e(imax - 1),b21d(imax),sigma21,dummy)

                 if (sigma11 <= sigma21) then
                    mu = sigma11
                    nu = sqrt(one - mu**2)
                    if (mu < thresh) then
                       mu = zero
                       nu = one
                    end if
                 else
                    nu = sigma21
                    mu = sqrt(1.0_sp - nu**2)
                    if (nu < thresh) then
                       mu = one
                       nu = zero
                    end if
                 end if
              end if
              ! rotate to produce bulges in b11 and b21
              if (mu <= nu) then
                 call la_slartgs(b11d(imin),b11e(imin),mu,work(iv1tcs + imin - 1),work(iv1tsn + &
                           imin - 1))
              else
                 call la_slartgs(b21d(imin),b21e(imin),nu,work(iv1tcs + imin - 1),work(iv1tsn + &
                           imin - 1))
              end if
              temp = work(iv1tcs + imin - 1)*b11d(imin) + work(iv1tsn + imin - 1)*b11e(imin)
              b11e(imin) = work(iv1tcs + imin - 1)*b11e(imin) - work(iv1tsn + imin - 1)*b11d(imin)
              b11d(imin) = temp
              b11bulge = work(iv1tsn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = work(iv1tcs + imin - 1)*b11d(imin + 1)
              temp = work(iv1tcs + imin - 1)*b21d(imin) + work(iv1tsn + imin - 1)*b21e(imin)
              b21e(imin) = work(iv1tcs + imin - 1)*b21e(imin) - work(iv1tsn + imin - 1)*b21d(imin)
              b21d(imin) = temp
              b21bulge = work(iv1tsn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = work(iv1tcs + imin - 1)*b21d(imin + 1)
              ! compute theta(imin)
              theta(imin) = atan2(sqrt(b21d(imin)**2 + b21bulge**2),sqrt(b11d(imin)**2 + &
                        b11bulge**2))
              ! chase the bulges in b11(imin+1,imin) and b21(imin+1,imin)
              if (b11d(imin)**2 + b11bulge**2 > thresh**2) then
                 call la_slartgp(b11bulge,b11d(imin),work(iu1sn + imin - 1),work(iu1cs + imin - 1), &
                            r)
              else if (mu <= nu) then
                 call la_slartgs(b11e(imin),b11d(imin + 1),mu,work(iu1cs + imin - 1),work( &
                           iu1sn + imin - 1))
              else
                 call la_slartgs(b12d(imin),b12e(imin),nu,work(iu1cs + imin - 1),work( &
                           iu1sn + imin - 1))
              end if
              if (b21d(imin)**2 + b21bulge**2 > thresh**2) then
                 call la_slartgp(b21bulge,b21d(imin),work(iu2sn + imin - 1),work(iu2cs + imin - 1), &
                            r)
              else if (nu < mu) then
                 call la_slartgs(b21e(imin),b21d(imin + 1),nu,work(iu2cs + imin - 1),work( &
                           iu2sn + imin - 1))
              else
                 call la_slartgs(b22d(imin),b22e(imin),mu,work(iu2cs + imin - 1),work(iu2sn + &
                           imin - 1))
              end if
              work(iu2cs + imin - 1) = -work(iu2cs + imin - 1)
              work(iu2sn + imin - 1) = -work(iu2sn + imin - 1)
              temp = work(iu1cs + imin - 1)*b11e(imin) + work(iu1sn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = work(iu1cs + imin - 1)*b11d(imin + 1) - work(iu1sn + imin - 1)*b11e(imin)

              b11e(imin) = temp
              if (imax > imin + 1) then
                 b11bulge = work(iu1sn + imin - 1)*b11e(imin + 1)
                 b11e(imin + 1) = work(iu1cs + imin - 1)*b11e(imin + 1)
              end if
              temp = work(iu1cs + imin - 1)*b12d(imin) + work(iu1sn + imin - 1)*b12e(imin)
              b12e(imin) = work(iu1cs + imin - 1)*b12e(imin) - work(iu1sn + imin - 1)*b12d(imin)
              b12d(imin) = temp
              b12bulge = work(iu1sn + imin - 1)*b12d(imin + 1)
              b12d(imin + 1) = work(iu1cs + imin - 1)*b12d(imin + 1)
              temp = work(iu2cs + imin - 1)*b21e(imin) + work(iu2sn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = work(iu2cs + imin - 1)*b21d(imin + 1) - work(iu2sn + imin - 1)*b21e(imin)

              b21e(imin) = temp
              if (imax > imin + 1) then
                 b21bulge = work(iu2sn + imin - 1)*b21e(imin + 1)
                 b21e(imin + 1) = work(iu2cs + imin - 1)*b21e(imin + 1)
              end if
              temp = work(iu2cs + imin - 1)*b22d(imin) + work(iu2sn + imin - 1)*b22e(imin)
              b22e(imin) = work(iu2cs + imin - 1)*b22e(imin) - work(iu2sn + imin - 1)*b22d(imin)
              b22d(imin) = temp
              b22bulge = work(iu2sn + imin - 1)*b22d(imin + 1)
              b22d(imin + 1) = work(iu2cs + imin - 1)*b22d(imin + 1)
              ! inner loop: chase bulges from b11(imin,imin+2),
              ! b12(imin,imin+1), b21(imin,imin+2), and b22(imin,imin+1) to
              ! bottom-right
              do i = imin + 1,imax - 1
                 ! compute phi(i-1)
                 x1 = sin(theta(i - 1))*b11e(i - 1) + cos(theta(i - 1))*b21e(i - 1)
                 x2 = sin(theta(i - 1))*b11bulge + cos(theta(i - 1))*b21bulge
                 y1 = sin(theta(i - 1))*b12d(i - 1) + cos(theta(i - 1))*b22d(i - 1)
                 y2 = sin(theta(i - 1))*b12bulge + cos(theta(i - 1))*b22bulge
                 phi(i - 1) = atan2(sqrt(x1**2 + x2**2),sqrt(y1**2 + y2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11e(i - 1)**2 + b11bulge**2 <= thresh**2
                 restart21 = b21e(i - 1)**2 + b21bulge**2 <= thresh**2
                 restart12 = b12d(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart22 = b22d(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i-1,i+1), b12(i-1,i),
                 ! b21(i-1,i+1), and b22(i-1,i). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart21) then
                    call la_slartgp(x2,x1,work(iv1tsn + i - 1),work(iv1tcs + i - 1),r)
                 else if (.not. restart11 .and. restart21) then
                    call la_slartgp(b11bulge,b11e(i - 1),work(iv1tsn + i - 1),work(iv1tcs + i - 1), &
                              r)
                 else if (restart11 .and. .not. restart21) then
                    call la_slartgp(b21bulge,b21e(i - 1),work(iv1tsn + i - 1),work(iv1tcs + i - 1), &
                              r)
                 else if (mu <= nu) then
                    call la_slartgs(b11d(i),b11e(i),mu,work(iv1tcs + i - 1),work(iv1tsn + i - 1))

                 else
                    call la_slartgs(b21d(i),b21e(i),nu,work(iv1tcs + i - 1),work(iv1tsn + i - 1))

                 end if
                 work(iv1tcs + i - 1) = -work(iv1tcs + i - 1)
                 work(iv1tsn + i - 1) = -work(iv1tsn + i - 1)
                 if (.not. restart12 .and. .not. restart22) then
                    call la_slartgp(y2,y1,work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - 1),r)

                 else if (.not. restart12 .and. restart22) then
                    call la_slartgp(b12bulge,b12d(i - 1),work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - &
                              1),r)
                 else if (restart12 .and. .not. restart22) then
                    call la_slartgp(b22bulge,b22d(i - 1),work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - &
                              1),r)
                 else if (nu < mu) then
                    call la_slartgs(b12e(i - 1),b12d(i),nu,work(iv2tcs + i - 1 - 1),work(iv2tsn + i - &
                              1 - 1))
                 else
                    call la_slartgs(b22e(i - 1),b22d(i),mu,work(iv2tcs + i - 1 - 1),work(iv2tsn + i - &
                              1 - 1))
                 end if
                 temp = work(iv1tcs + i - 1)*b11d(i) + work(iv1tsn + i - 1)*b11e(i)
                 b11e(i) = work(iv1tcs + i - 1)*b11e(i) - work(iv1tsn + i - 1)*b11d(i)
                 b11d(i) = temp
                 b11bulge = work(iv1tsn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = work(iv1tcs + i - 1)*b11d(i + 1)
                 temp = work(iv1tcs + i - 1)*b21d(i) + work(iv1tsn + i - 1)*b21e(i)
                 b21e(i) = work(iv1tcs + i - 1)*b21e(i) - work(iv1tsn + i - 1)*b21d(i)
                 b21d(i) = temp
                 b21bulge = work(iv1tsn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = work(iv1tcs + i - 1)*b21d(i + 1)
                 temp = work(iv2tcs + i - 1 - 1)*b12e(i - 1) + work(iv2tsn + i - 1 - 1)*b12d(i)
                 b12d(i) = work(iv2tcs + i - 1 - 1)*b12d(i) - work(iv2tsn + i - 1 - 1)*b12e(i - 1)
                 b12e(i - 1) = temp
                 b12bulge = work(iv2tsn + i - 1 - 1)*b12e(i)
                 b12e(i) = work(iv2tcs + i - 1 - 1)*b12e(i)
                 temp = work(iv2tcs + i - 1 - 1)*b22e(i - 1) + work(iv2tsn + i - 1 - 1)*b22d(i)
                 b22d(i) = work(iv2tcs + i - 1 - 1)*b22d(i) - work(iv2tsn + i - 1 - 1)*b22e(i - 1)
                 b22e(i - 1) = temp
                 b22bulge = work(iv2tsn + i - 1 - 1)*b22e(i)
                 b22e(i) = work(iv2tcs + i - 1 - 1)*b22e(i)
                 ! compute theta(i)
                 x1 = cos(phi(i - 1))*b11d(i) + sin(phi(i - 1))*b12e(i - 1)
                 x2 = cos(phi(i - 1))*b11bulge + sin(phi(i - 1))*b12bulge
                 y1 = cos(phi(i - 1))*b21d(i) + sin(phi(i - 1))*b22e(i - 1)
                 y2 = cos(phi(i - 1))*b21bulge + sin(phi(i - 1))*b22bulge
                 theta(i) = atan2(sqrt(y1**2 + y2**2),sqrt(x1**2 + x2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11d(i)**2 + b11bulge**2 <= thresh**2
                 restart12 = b12e(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart21 = b21d(i)**2 + b21bulge**2 <= thresh**2
                 restart22 = b22e(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i+1,i), b12(i+1,i-1),
                 ! b21(i+1,i), and b22(i+1,i-1). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart12) then
                    call la_slartgp(x2,x1,work(iu1sn + i - 1),work(iu1cs + i - 1),r)
                 else if (.not. restart11 .and. restart12) then
                    call la_slartgp(b11bulge,b11d(i),work(iu1sn + i - 1),work(iu1cs + i - 1),r)

                 else if (restart11 .and. .not. restart12) then
                    call la_slartgp(b12bulge,b12e(i - 1),work(iu1sn + i - 1),work(iu1cs + i - 1),r)

                 else if (mu <= nu) then
                    call la_slartgs(b11e(i),b11d(i + 1),mu,work(iu1cs + i - 1),work(iu1sn + i - 1))

                 else
                    call la_slartgs(b12d(i),b12e(i),nu,work(iu1cs + i - 1),work(iu1sn + i - 1))

                 end if
                 if (.not. restart21 .and. .not. restart22) then
                    call la_slartgp(y2,y1,work(iu2sn + i - 1),work(iu2cs + i - 1),r)
                 else if (.not. restart21 .and. restart22) then
                    call la_slartgp(b21bulge,b21d(i),work(iu2sn + i - 1),work(iu2cs + i - 1),r)

                 else if (restart21 .and. .not. restart22) then
                    call la_slartgp(b22bulge,b22e(i - 1),work(iu2sn + i - 1),work(iu2cs + i - 1),r)

                 else if (nu < mu) then
                    call la_slartgs(b21e(i),b21e(i + 1),nu,work(iu2cs + i - 1),work(iu2sn + i - 1))

                 else
                    call la_slartgs(b22d(i),b22e(i),mu,work(iu2cs + i - 1),work(iu2sn + i - 1))

                 end if
                 work(iu2cs + i - 1) = -work(iu2cs + i - 1)
                 work(iu2sn + i - 1) = -work(iu2sn + i - 1)
                 temp = work(iu1cs + i - 1)*b11e(i) + work(iu1sn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = work(iu1cs + i - 1)*b11d(i + 1) - work(iu1sn + i - 1)*b11e(i)
                 b11e(i) = temp
                 if (i < imax - 1) then
                    b11bulge = work(iu1sn + i - 1)*b11e(i + 1)
                    b11e(i + 1) = work(iu1cs + i - 1)*b11e(i + 1)
                 end if
                 temp = work(iu2cs + i - 1)*b21e(i) + work(iu2sn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = work(iu2cs + i - 1)*b21d(i + 1) - work(iu2sn + i - 1)*b21e(i)
                 b21e(i) = temp
                 if (i < imax - 1) then
                    b21bulge = work(iu2sn + i - 1)*b21e(i + 1)
                    b21e(i + 1) = work(iu2cs + i - 1)*b21e(i + 1)
                 end if
                 temp = work(iu1cs + i - 1)*b12d(i) + work(iu1sn + i - 1)*b12e(i)
                 b12e(i) = work(iu1cs + i - 1)*b12e(i) - work(iu1sn + i - 1)*b12d(i)
                 b12d(i) = temp
                 b12bulge = work(iu1sn + i - 1)*b12d(i + 1)
                 b12d(i + 1) = work(iu1cs + i - 1)*b12d(i + 1)
                 temp = work(iu2cs + i - 1)*b22d(i) + work(iu2sn + i - 1)*b22e(i)
                 b22e(i) = work(iu2cs + i - 1)*b22e(i) - work(iu2sn + i - 1)*b22d(i)
                 b22d(i) = temp
                 b22bulge = work(iu2sn + i - 1)*b22d(i + 1)
                 b22d(i + 1) = work(iu2cs + i - 1)*b22d(i + 1)
              end do
              ! compute phi(imax-1)
              x1 = sin(theta(imax - 1))*b11e(imax - 1) + cos(theta(imax - 1))*b21e(imax - 1)
              y1 = sin(theta(imax - 1))*b12d(imax - 1) + cos(theta(imax - 1))*b22d(imax - 1)
              y2 = sin(theta(imax - 1))*b12bulge + cos(theta(imax - 1))*b22bulge
              phi(imax - 1) = atan2(abs(x1),sqrt(y1**2 + y2**2))
              ! chase bulges from b12(imax-1,imax) and b22(imax-1,imax)
              restart12 = b12d(imax - 1)**2 + b12bulge**2 <= thresh**2
              restart22 = b22d(imax - 1)**2 + b22bulge**2 <= thresh**2
              if (.not. restart12 .and. .not. restart22) then
                 call la_slartgp(y2,y1,work(iv2tsn + imax - 1 - 1),work(iv2tcs + imax - 1 - 1),r)

              else if (.not. restart12 .and. restart22) then
                 call la_slartgp(b12bulge,b12d(imax - 1),work(iv2tsn + imax - 1 - 1),work(iv2tcs + &
                           imax - 1 - 1),r)
              else if (restart12 .and. .not. restart22) then
                 call la_slartgp(b22bulge,b22d(imax - 1),work(iv2tsn + imax - 1 - 1),work(iv2tcs + &
                           imax - 1 - 1),r)
              else if (nu < mu) then
                 call la_slartgs(b12e(imax - 1),b12d(imax),nu,work(iv2tcs + imax - 1 - 1),work( &
                           iv2tsn + imax - 1 - 1))
              else
                 call la_slartgs(b22e(imax - 1),b22d(imax),mu,work(iv2tcs + imax - 1 - 1),work( &
                           iv2tsn + imax - 1 - 1))
              end if
              temp = work(iv2tcs + imax - 1 - 1)*b12e(imax - 1) + work(iv2tsn + imax - 1 - 1)*b12d(imax)
              b12d(imax) = work(iv2tcs + imax - 1 - 1)*b12d(imax) - work(iv2tsn + imax - 1 - 1)*b12e(imax - 1)

              b12e(imax - 1) = temp
              temp = work(iv2tcs + imax - 1 - 1)*b22e(imax - 1) + work(iv2tsn + imax - 1 - 1)*b22d(imax)
              b22d(imax) = work(iv2tcs + imax - 1 - 1)*b22d(imax) - work(iv2tsn + imax - 1 - 1)*b22e(imax - 1)

              b22e(imax - 1) = temp
              ! update singular vectors
              if (wantu1) then
                 if (colmajor) then
                    call la_slasr('R','V','F',p,imax - imin + 1,work(iu1cs + imin - 1),work( &
                              iu1sn + imin - 1),u1(1,imin),ldu1)
                 else
                    call la_slasr('L','V','F',imax - imin + 1,p,work(iu1cs + imin - 1),work( &
                              iu1sn + imin - 1),u1(imin,1),ldu1)
                 end if
              end if
              if (wantu2) then
                 if (colmajor) then
                    call la_slasr('R','V','F',m - p,imax - imin + 1,work(iu2cs + imin - 1),work( &
                              iu2sn + imin - 1),u2(1,imin),ldu2)
                 else
                    call la_slasr('L','V','F',imax - imin + 1,m - p,work(iu2cs + imin - 1),work( &
                              iu2sn + imin - 1),u2(imin,1),ldu2)
                 end if
              end if
              if (wantv1t) then
                 if (colmajor) then
                    call la_slasr('L','V','F',imax - imin + 1,q,work(iv1tcs + imin - 1),work( &
                              iv1tsn + imin - 1),v1t(imin,1),ldv1t)
                 else
                    call la_slasr('R','V','F',q,imax - imin + 1,work(iv1tcs + imin - 1),work( &
                              iv1tsn + imin - 1),v1t(1,imin),ldv1t)
                 end if
              end if
              if (wantv2t) then
                 if (colmajor) then
                    call la_slasr('L','V','F',imax - imin + 1,m - q,work(iv2tcs + imin - 1),work( &
                              iv2tsn + imin - 1),v2t(imin,1),ldv2t)
                 else
                    call la_slasr('R','V','F',m - q,imax - imin + 1,work(iv2tcs + imin - 1),work( &
                              iv2tsn + imin - 1),v2t(1,imin),ldv2t)
                 end if
              end if
              ! fix signs on b11(imax-1,imax) and b21(imax-1,imax)
              if (b11e(imax - 1) + b21e(imax - 1) > 0) then
                 b11d(imax) = -b11d(imax)
                 b21d(imax) = -b21d(imax)
                 if (wantv1t) then
                    if (colmajor) then
                       call la_sscal(q,negone,v1t(imax,1),ldv1t)
                    else
                       call la_sscal(q,negone,v1t(1,imax),1)
                    end if
                 end if
              end if
              ! compute theta(imax)
              x1 = cos(phi(imax - 1))*b11d(imax) + sin(phi(imax - 1))*b12e(imax - 1)
              y1 = cos(phi(imax - 1))*b21d(imax) + sin(phi(imax - 1))*b22e(imax - 1)
              theta(imax) = atan2(abs(y1),abs(x1))
              ! fix signs on b11(imax,imax), b12(imax,imax-1), b21(imax,imax),
              ! and b22(imax,imax-1)
              if (b11d(imax) + b12e(imax - 1) < 0) then
                 b12d(imax) = -b12d(imax)
                 if (wantu1) then
                    if (colmajor) then
                       call la_sscal(p,negone,u1(1,imax),1)
                    else
                       call la_sscal(p,negone,u1(imax,1),ldu1)
                    end if
                 end if
              end if
              if (b21d(imax) + b22e(imax - 1) > 0) then
                 b22d(imax) = -b22d(imax)
                 if (wantu2) then
                    if (colmajor) then
                       call la_sscal(m - p,negone,u2(1,imax),1)
                    else
                       call la_sscal(m - p,negone,u2(imax,1),ldu2)
                    end if
                 end if
              end if
              ! fix signs on b12(imax,imax) and b22(imax,imax)
              if (b12d(imax) + b22d(imax) < 0) then
                 if (wantv2t) then
                    if (colmajor) then
                       call la_sscal(m - q,negone,v2t(imax,1),ldv2t)
                    else
                       call la_sscal(m - q,negone,v2t(1,imax),1)
                    end if
                 end if
              end if
              ! test for negligible sines or cosines
              do i = imin,imax
                 if (theta(i) < thresh) then
                    theta(i) = zero
                 else if (theta(i) > piover2 - thresh) then
                    theta(i) = piover2
                 end if
              end do
              do i = imin,imax - 1
                 if (phi(i) < thresh) then
                    phi(i) = zero
                 else if (phi(i) > piover2 - thresh) then
                    phi(i) = piover2
                 end if
              end do
              ! deflate
              if (imax > 1) then
                 do while (phi(imax - 1) == zero)
                    imax = imax - 1
                    if (imax <= 1) exit
                 end do
              end if
              if (imin > imax - 1) imin = imax - 1
              if (imin > 1) then
                 do while (phi(imin - 1) /= zero)
                     imin = imin - 1
                     if (imin <= 1) exit
                 end do
              end if
              ! repeat main iteration loop
           end do
           ! postprocessing: order theta from least to greatest
           do i = 1,q
              mini = i
              thetamin = theta(i)
              do j = i + 1,q
                 if (theta(j) < thetamin) then
                    mini = j
                    thetamin = theta(j)
                 end if
              end do
              if (mini /= i) then
                 theta(mini) = theta(i)
                 theta(i) = thetamin
                 if (colmajor) then
                    if (wantu1) call la_sswap(p,u1(1,i),1,u1(1,mini),1)
                    if (wantu2) call la_sswap(m - p,u2(1,i),1,u2(1,mini),1)
                    if (wantv1t) call la_sswap(q,v1t(i,1),ldv1t,v1t(mini,1),ldv1t)

                    if (wantv2t) call la_sswap(m - q,v2t(i,1),ldv2t,v2t(mini,1),ldv2t)

                 else
                    if (wantu1) call la_sswap(p,u1(i,1),ldu1,u1(mini,1),ldu1)
                    if (wantu2) call la_sswap(m - p,u2(i,1),ldu2,u2(mini,1),ldu2)
                    if (wantv1t) call la_sswap(q,v1t(1,i),1,v1t(1,mini),1)
                    if (wantv2t) call la_sswap(m - q,v2t(1,i),1,v2t(1,mini),1)
                 end if
              end if
           end do
           return
     end subroutine la_sbbcsd
     !> DBBCSD: computes the CS decomposition of an orthogonal matrix in
     !> bidiagonal-block form,
     !> [ B11 | B12 0  0 ]
     !> [  0  |  0 -I  0 ]
     !> X = [----------------]
     !> [ B21 | B22 0  0 ]
     !> [  0  |  0  0  I ]
     !> [  C | -S  0  0 ]
     !> [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T
     !> = [---------] [---------------] [---------]   .
     !> [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
     !> [  0 |  0  0  I ]
     !> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
     !> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
     !> transposed and/or permuted. This can be done in constant time using
     !> the TRANS and SIGNS options. See DORCSD for details.)
     !> The bidiagonal matrices B11, B12, B21, and B22 are represented
     !> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
     !> The orthogonal matrices U1, U2, V1T, and V2T are input/output.
     !> The input matrices are pre- or post-multiplied by the appropriate
     !> singular vector matrices.

     pure subroutine la_dbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,phi,u1, &
     ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,b11d,b11e,b12d,b12e,b21d,b21e,b22d,b22e,work, &
               lwork,info)
        use la_constants_dp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,lwork,m,p,q
           ! Array Arguments
           real(dp),intent(out) :: b11d(*),b11e(*),b12d(*),b12e(*),b21d(*),b21e(*),b22d(*), &
                      b22e(*),work(*)
           real(dp),intent(inout) :: phi(*),theta(*)
           real(dp),intent(inout) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*)
        ! ===================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 6
           real(dp),parameter :: hundred = 100.0_dp
           real(dp),parameter :: meighth = -0.125_dp
           real(dp),parameter :: piover2 = 1.57079632679489661923132169163975144210_dp

           ! Local Scalars
           logical(lk) :: colmajor,lquery,restart11,restart12,restart21,restart22,wantu1, &
                     wantu2,wantv1t,wantv2t
           integer(ilp) :: i,imin,imax,iter,iu1cs,iu1sn,iu2cs,iu2sn,iv1tcs,iv1tsn, &
                     iv2tcs,iv2tsn,j,lworkmin,lworkopt,maxit,mini
           real(dp) :: b11bulge,b12bulge,b21bulge,b22bulge,dummy,eps,mu,nu,r,sigma11, &
                     sigma21,temp,thetamax,thetamin,thresh,tol,tolmul,unfl,x1,x2,y1,y2
           ! Intrinsic Functions
           intrinsic :: abs,atan2,cos,max,min,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           if (m < 0) then
              info = -6
           else if (p < 0 .or. p > m) then
              info = -7
           else if (q < 0 .or. q > m) then
              info = -8
           else if (q > p .or. q > m - p .or. q > m - q) then
              info = -8
           else if (wantu1 .and. ldu1 < p) then
              info = -12
           else if (wantu2 .and. ldu2 < m - p) then
              info = -14
           else if (wantv1t .and. ldv1t < q) then
              info = -16
           else if (wantv2t .and. ldv2t < m - q) then
              info = -18
           end if
           ! quick return if q = 0
           if (info == 0 .and. q == 0) then
              lworkmin = 1
              work(1) = lworkmin
              return
           end if
           ! compute workspace
           if (info == 0) then
              iu1cs = 1
              iu1sn = iu1cs + q
              iu2cs = iu1sn + q
              iu2sn = iu2cs + q
              iv1tcs = iu2sn + q
              iv1tsn = iv1tcs + q
              iv2tcs = iv1tsn + q
              iv2tsn = iv2tcs + q
              lworkopt = iv2tsn + q - 1
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -28
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DBBCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! get machine constants
           eps = la_dlamch('EPSILON')
           unfl = la_dlamch('SAFE MINIMUM')
           tolmul = max(ten,min(hundred,eps**meighth))
           tol = tolmul*eps
           thresh = max(tol,maxitr*q*q*unfl)
           ! test for negligible sines or cosines
           do i = 1,q
              if (theta(i) < thresh) then
                 theta(i) = zero
              else if (theta(i) > piover2 - thresh) then
                 theta(i) = piover2
              end if
           end do
           do i = 1,q - 1
              if (phi(i) < thresh) then
                 phi(i) = zero
              else if (phi(i) > piover2 - thresh) then
                 phi(i) = piover2
              end if
           end do
           ! initial deflation
           imax = q
           do while (imax > 1)
              if (phi(imax - 1) /= zero) then
                 exit
              end if
              imax = imax - 1
           end do
           imin = imax - 1
           if (imin > 1) then
              do while (phi(imin - 1) /= zero)
                 imin = imin - 1
                 if (imin <= 1) exit
              end do
           end if
           ! initialize iteration counter
           maxit = maxitr*q*q
           iter = 0
           ! begin main iteration loop
           do while (imax > 1)
              ! compute the matrix entries
              b11d(imin) = cos(theta(imin))
              b21d(imin) = -sin(theta(imin))
              do i = imin,imax - 1
                 b11e(i) = -sin(theta(i))*sin(phi(i))
                 b11d(i + 1) = cos(theta(i + 1))*cos(phi(i))
                 b12d(i) = sin(theta(i))*cos(phi(i))
                 b12e(i) = cos(theta(i + 1))*sin(phi(i))
                 b21e(i) = -cos(theta(i))*sin(phi(i))
                 b21d(i + 1) = -sin(theta(i + 1))*cos(phi(i))
                 b22d(i) = cos(theta(i))*cos(phi(i))
                 b22e(i) = -sin(theta(i + 1))*sin(phi(i))
              end do
              b12d(imax) = sin(theta(imax))
              b22d(imax) = cos(theta(imax))
              ! abort if not converging; otherwise, increment iter
              if (iter > maxit) then
                 info = 0
                 do i = 1,q
                    if (phi(i) /= zero) info = info + 1
                 end do
                 return
              end if
              iter = iter + imax - imin
              ! compute shifts
              thetamax = theta(imin)
              thetamin = theta(imin)
              do i = imin + 1,imax
                 if (theta(i) > thetamax) thetamax = theta(i)
                 if (theta(i) < thetamin) thetamin = theta(i)
              end do
              if (thetamax > piover2 - thresh) then
                 ! zero on diagonals of b11 and b22; induce deflation with a
                 ! zero shift
                 mu = zero
                 nu = one
              else if (thetamin < thresh) then
                 ! zero on diagonals of b12 and b22; induce deflation with a
                 ! zero shift
                 mu = one
                 nu = zero
              else
                 ! compute shifts for b11 and b21 and use the lesser
                 call la_dlas2(b11d(imax - 1),b11e(imax - 1),b11d(imax),sigma11,dummy)

                 call la_dlas2(b21d(imax - 1),b21e(imax - 1),b21d(imax),sigma21,dummy)

                 if (sigma11 <= sigma21) then
                    mu = sigma11
                    nu = sqrt(one - mu**2)
                    if (mu < thresh) then
                       mu = zero
                       nu = one
                    end if
                 else
                    nu = sigma21
                    mu = sqrt(1.0_dp - nu**2)
                    if (nu < thresh) then
                       mu = one
                       nu = zero
                    end if
                 end if
              end if
              ! rotate to produce bulges in b11 and b21
              if (mu <= nu) then
                 call la_dlartgs(b11d(imin),b11e(imin),mu,work(iv1tcs + imin - 1),work(iv1tsn + &
                           imin - 1))
              else
                 call la_dlartgs(b21d(imin),b21e(imin),nu,work(iv1tcs + imin - 1),work(iv1tsn + &
                           imin - 1))
              end if
              temp = work(iv1tcs + imin - 1)*b11d(imin) + work(iv1tsn + imin - 1)*b11e(imin)
              b11e(imin) = work(iv1tcs + imin - 1)*b11e(imin) - work(iv1tsn + imin - 1)*b11d(imin)
              b11d(imin) = temp
              b11bulge = work(iv1tsn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = work(iv1tcs + imin - 1)*b11d(imin + 1)
              temp = work(iv1tcs + imin - 1)*b21d(imin) + work(iv1tsn + imin - 1)*b21e(imin)
              b21e(imin) = work(iv1tcs + imin - 1)*b21e(imin) - work(iv1tsn + imin - 1)*b21d(imin)
              b21d(imin) = temp
              b21bulge = work(iv1tsn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = work(iv1tcs + imin - 1)*b21d(imin + 1)
              ! compute theta(imin)
              theta(imin) = atan2(sqrt(b21d(imin)**2 + b21bulge**2),sqrt(b11d(imin)**2 + &
                        b11bulge**2))
              ! chase the bulges in b11(imin+1,imin) and b21(imin+1,imin)
              if (b11d(imin)**2 + b11bulge**2 > thresh**2) then
                 call la_dlartgp(b11bulge,b11d(imin),work(iu1sn + imin - 1),work(iu1cs + imin - 1), &
                            r)
              else if (mu <= nu) then
                 call la_dlartgs(b11e(imin),b11d(imin + 1),mu,work(iu1cs + imin - 1),work( &
                           iu1sn + imin - 1))
              else
                 call la_dlartgs(b12d(imin),b12e(imin),nu,work(iu1cs + imin - 1),work( &
                           iu1sn + imin - 1))
              end if
              if (b21d(imin)**2 + b21bulge**2 > thresh**2) then
                 call la_dlartgp(b21bulge,b21d(imin),work(iu2sn + imin - 1),work(iu2cs + imin - 1), &
                            r)
              else if (nu < mu) then
                 call la_dlartgs(b21e(imin),b21d(imin + 1),nu,work(iu2cs + imin - 1),work( &
                           iu2sn + imin - 1))
              else
                 call la_dlartgs(b22d(imin),b22e(imin),mu,work(iu2cs + imin - 1),work(iu2sn + &
                           imin - 1))
              end if
              work(iu2cs + imin - 1) = -work(iu2cs + imin - 1)
              work(iu2sn + imin - 1) = -work(iu2sn + imin - 1)
              temp = work(iu1cs + imin - 1)*b11e(imin) + work(iu1sn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = work(iu1cs + imin - 1)*b11d(imin + 1) - work(iu1sn + imin - 1)*b11e(imin)

              b11e(imin) = temp
              if (imax > imin + 1) then
                 b11bulge = work(iu1sn + imin - 1)*b11e(imin + 1)
                 b11e(imin + 1) = work(iu1cs + imin - 1)*b11e(imin + 1)
              end if
              temp = work(iu1cs + imin - 1)*b12d(imin) + work(iu1sn + imin - 1)*b12e(imin)
              b12e(imin) = work(iu1cs + imin - 1)*b12e(imin) - work(iu1sn + imin - 1)*b12d(imin)
              b12d(imin) = temp
              b12bulge = work(iu1sn + imin - 1)*b12d(imin + 1)
              b12d(imin + 1) = work(iu1cs + imin - 1)*b12d(imin + 1)
              temp = work(iu2cs + imin - 1)*b21e(imin) + work(iu2sn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = work(iu2cs + imin - 1)*b21d(imin + 1) - work(iu2sn + imin - 1)*b21e(imin)

              b21e(imin) = temp
              if (imax > imin + 1) then
                 b21bulge = work(iu2sn + imin - 1)*b21e(imin + 1)
                 b21e(imin + 1) = work(iu2cs + imin - 1)*b21e(imin + 1)
              end if
              temp = work(iu2cs + imin - 1)*b22d(imin) + work(iu2sn + imin - 1)*b22e(imin)
              b22e(imin) = work(iu2cs + imin - 1)*b22e(imin) - work(iu2sn + imin - 1)*b22d(imin)
              b22d(imin) = temp
              b22bulge = work(iu2sn + imin - 1)*b22d(imin + 1)
              b22d(imin + 1) = work(iu2cs + imin - 1)*b22d(imin + 1)
              ! inner loop: chase bulges from b11(imin,imin+2),
              ! b12(imin,imin+1), b21(imin,imin+2), and b22(imin,imin+1) to
              ! bottom-right
              do i = imin + 1,imax - 1
                 ! compute phi(i-1)
                 x1 = sin(theta(i - 1))*b11e(i - 1) + cos(theta(i - 1))*b21e(i - 1)
                 x2 = sin(theta(i - 1))*b11bulge + cos(theta(i - 1))*b21bulge
                 y1 = sin(theta(i - 1))*b12d(i - 1) + cos(theta(i - 1))*b22d(i - 1)
                 y2 = sin(theta(i - 1))*b12bulge + cos(theta(i - 1))*b22bulge
                 phi(i - 1) = atan2(sqrt(x1**2 + x2**2),sqrt(y1**2 + y2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11e(i - 1)**2 + b11bulge**2 <= thresh**2
                 restart21 = b21e(i - 1)**2 + b21bulge**2 <= thresh**2
                 restart12 = b12d(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart22 = b22d(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i-1,i+1), b12(i-1,i),
                 ! b21(i-1,i+1), and b22(i-1,i). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart21) then
                    call la_dlartgp(x2,x1,work(iv1tsn + i - 1),work(iv1tcs + i - 1),r)
                 else if (.not. restart11 .and. restart21) then
                    call la_dlartgp(b11bulge,b11e(i - 1),work(iv1tsn + i - 1),work(iv1tcs + i - 1), &
                              r)
                 else if (restart11 .and. .not. restart21) then
                    call la_dlartgp(b21bulge,b21e(i - 1),work(iv1tsn + i - 1),work(iv1tcs + i - 1), &
                              r)
                 else if (mu <= nu) then
                    call la_dlartgs(b11d(i),b11e(i),mu,work(iv1tcs + i - 1),work(iv1tsn + i - 1))

                 else
                    call la_dlartgs(b21d(i),b21e(i),nu,work(iv1tcs + i - 1),work(iv1tsn + i - 1))

                 end if
                 work(iv1tcs + i - 1) = -work(iv1tcs + i - 1)
                 work(iv1tsn + i - 1) = -work(iv1tsn + i - 1)
                 if (.not. restart12 .and. .not. restart22) then
                    call la_dlartgp(y2,y1,work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - 1),r)

                 else if (.not. restart12 .and. restart22) then
                    call la_dlartgp(b12bulge,b12d(i - 1),work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - &
                              1),r)
                 else if (restart12 .and. .not. restart22) then
                    call la_dlartgp(b22bulge,b22d(i - 1),work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - &
                              1),r)
                 else if (nu < mu) then
                    call la_dlartgs(b12e(i - 1),b12d(i),nu,work(iv2tcs + i - 1 - 1),work(iv2tsn + i - &
                              1 - 1))
                 else
                    call la_dlartgs(b22e(i - 1),b22d(i),mu,work(iv2tcs + i - 1 - 1),work(iv2tsn + i - &
                              1 - 1))
                 end if
                 temp = work(iv1tcs + i - 1)*b11d(i) + work(iv1tsn + i - 1)*b11e(i)
                 b11e(i) = work(iv1tcs + i - 1)*b11e(i) - work(iv1tsn + i - 1)*b11d(i)
                 b11d(i) = temp
                 b11bulge = work(iv1tsn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = work(iv1tcs + i - 1)*b11d(i + 1)
                 temp = work(iv1tcs + i - 1)*b21d(i) + work(iv1tsn + i - 1)*b21e(i)
                 b21e(i) = work(iv1tcs + i - 1)*b21e(i) - work(iv1tsn + i - 1)*b21d(i)
                 b21d(i) = temp
                 b21bulge = work(iv1tsn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = work(iv1tcs + i - 1)*b21d(i + 1)
                 temp = work(iv2tcs + i - 1 - 1)*b12e(i - 1) + work(iv2tsn + i - 1 - 1)*b12d(i)
                 b12d(i) = work(iv2tcs + i - 1 - 1)*b12d(i) - work(iv2tsn + i - 1 - 1)*b12e(i - 1)
                 b12e(i - 1) = temp
                 b12bulge = work(iv2tsn + i - 1 - 1)*b12e(i)
                 b12e(i) = work(iv2tcs + i - 1 - 1)*b12e(i)
                 temp = work(iv2tcs + i - 1 - 1)*b22e(i - 1) + work(iv2tsn + i - 1 - 1)*b22d(i)
                 b22d(i) = work(iv2tcs + i - 1 - 1)*b22d(i) - work(iv2tsn + i - 1 - 1)*b22e(i - 1)
                 b22e(i - 1) = temp
                 b22bulge = work(iv2tsn + i - 1 - 1)*b22e(i)
                 b22e(i) = work(iv2tcs + i - 1 - 1)*b22e(i)
                 ! compute theta(i)
                 x1 = cos(phi(i - 1))*b11d(i) + sin(phi(i - 1))*b12e(i - 1)
                 x2 = cos(phi(i - 1))*b11bulge + sin(phi(i - 1))*b12bulge
                 y1 = cos(phi(i - 1))*b21d(i) + sin(phi(i - 1))*b22e(i - 1)
                 y2 = cos(phi(i - 1))*b21bulge + sin(phi(i - 1))*b22bulge
                 theta(i) = atan2(sqrt(y1**2 + y2**2),sqrt(x1**2 + x2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11d(i)**2 + b11bulge**2 <= thresh**2
                 restart12 = b12e(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart21 = b21d(i)**2 + b21bulge**2 <= thresh**2
                 restart22 = b22e(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i+1,i), b12(i+1,i-1),
                 ! b21(i+1,i), and b22(i+1,i-1). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart12) then
                    call la_dlartgp(x2,x1,work(iu1sn + i - 1),work(iu1cs + i - 1),r)
                 else if (.not. restart11 .and. restart12) then
                    call la_dlartgp(b11bulge,b11d(i),work(iu1sn + i - 1),work(iu1cs + i - 1),r)

                 else if (restart11 .and. .not. restart12) then
                    call la_dlartgp(b12bulge,b12e(i - 1),work(iu1sn + i - 1),work(iu1cs + i - 1),r)

                 else if (mu <= nu) then
                    call la_dlartgs(b11e(i),b11d(i + 1),mu,work(iu1cs + i - 1),work(iu1sn + i - 1))

                 else
                    call la_dlartgs(b12d(i),b12e(i),nu,work(iu1cs + i - 1),work(iu1sn + i - 1))

                 end if
                 if (.not. restart21 .and. .not. restart22) then
                    call la_dlartgp(y2,y1,work(iu2sn + i - 1),work(iu2cs + i - 1),r)
                 else if (.not. restart21 .and. restart22) then
                    call la_dlartgp(b21bulge,b21d(i),work(iu2sn + i - 1),work(iu2cs + i - 1),r)

                 else if (restart21 .and. .not. restart22) then
                    call la_dlartgp(b22bulge,b22e(i - 1),work(iu2sn + i - 1),work(iu2cs + i - 1),r)

                 else if (nu < mu) then
                    call la_dlartgs(b21e(i),b21e(i + 1),nu,work(iu2cs + i - 1),work(iu2sn + i - 1))

                 else
                    call la_dlartgs(b22d(i),b22e(i),mu,work(iu2cs + i - 1),work(iu2sn + i - 1))

                 end if
                 work(iu2cs + i - 1) = -work(iu2cs + i - 1)
                 work(iu2sn + i - 1) = -work(iu2sn + i - 1)
                 temp = work(iu1cs + i - 1)*b11e(i) + work(iu1sn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = work(iu1cs + i - 1)*b11d(i + 1) - work(iu1sn + i - 1)*b11e(i)
                 b11e(i) = temp
                 if (i < imax - 1) then
                    b11bulge = work(iu1sn + i - 1)*b11e(i + 1)
                    b11e(i + 1) = work(iu1cs + i - 1)*b11e(i + 1)
                 end if
                 temp = work(iu2cs + i - 1)*b21e(i) + work(iu2sn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = work(iu2cs + i - 1)*b21d(i + 1) - work(iu2sn + i - 1)*b21e(i)
                 b21e(i) = temp
                 if (i < imax - 1) then
                    b21bulge = work(iu2sn + i - 1)*b21e(i + 1)
                    b21e(i + 1) = work(iu2cs + i - 1)*b21e(i + 1)
                 end if
                 temp = work(iu1cs + i - 1)*b12d(i) + work(iu1sn + i - 1)*b12e(i)
                 b12e(i) = work(iu1cs + i - 1)*b12e(i) - work(iu1sn + i - 1)*b12d(i)
                 b12d(i) = temp
                 b12bulge = work(iu1sn + i - 1)*b12d(i + 1)
                 b12d(i + 1) = work(iu1cs + i - 1)*b12d(i + 1)
                 temp = work(iu2cs + i - 1)*b22d(i) + work(iu2sn + i - 1)*b22e(i)
                 b22e(i) = work(iu2cs + i - 1)*b22e(i) - work(iu2sn + i - 1)*b22d(i)
                 b22d(i) = temp
                 b22bulge = work(iu2sn + i - 1)*b22d(i + 1)
                 b22d(i + 1) = work(iu2cs + i - 1)*b22d(i + 1)
              end do
              ! compute phi(imax-1)
              x1 = sin(theta(imax - 1))*b11e(imax - 1) + cos(theta(imax - 1))*b21e(imax - 1)
              y1 = sin(theta(imax - 1))*b12d(imax - 1) + cos(theta(imax - 1))*b22d(imax - 1)
              y2 = sin(theta(imax - 1))*b12bulge + cos(theta(imax - 1))*b22bulge
              phi(imax - 1) = atan2(abs(x1),sqrt(y1**2 + y2**2))
              ! chase bulges from b12(imax-1,imax) and b22(imax-1,imax)
              restart12 = b12d(imax - 1)**2 + b12bulge**2 <= thresh**2
              restart22 = b22d(imax - 1)**2 + b22bulge**2 <= thresh**2
              if (.not. restart12 .and. .not. restart22) then
                 call la_dlartgp(y2,y1,work(iv2tsn + imax - 1 - 1),work(iv2tcs + imax - 1 - 1),r)

              else if (.not. restart12 .and. restart22) then
                 call la_dlartgp(b12bulge,b12d(imax - 1),work(iv2tsn + imax - 1 - 1),work(iv2tcs + &
                           imax - 1 - 1),r)
              else if (restart12 .and. .not. restart22) then
                 call la_dlartgp(b22bulge,b22d(imax - 1),work(iv2tsn + imax - 1 - 1),work(iv2tcs + &
                           imax - 1 - 1),r)
              else if (nu < mu) then
                 call la_dlartgs(b12e(imax - 1),b12d(imax),nu,work(iv2tcs + imax - 1 - 1),work( &
                           iv2tsn + imax - 1 - 1))
              else
                 call la_dlartgs(b22e(imax - 1),b22d(imax),mu,work(iv2tcs + imax - 1 - 1),work( &
                           iv2tsn + imax - 1 - 1))
              end if
              temp = work(iv2tcs + imax - 1 - 1)*b12e(imax - 1) + work(iv2tsn + imax - 1 - 1)*b12d(imax)
              b12d(imax) = work(iv2tcs + imax - 1 - 1)*b12d(imax) - work(iv2tsn + imax - 1 - 1)*b12e(imax - 1)

              b12e(imax - 1) = temp
              temp = work(iv2tcs + imax - 1 - 1)*b22e(imax - 1) + work(iv2tsn + imax - 1 - 1)*b22d(imax)
              b22d(imax) = work(iv2tcs + imax - 1 - 1)*b22d(imax) - work(iv2tsn + imax - 1 - 1)*b22e(imax - 1)

              b22e(imax - 1) = temp
              ! update singular vectors
              if (wantu1) then
                 if (colmajor) then
                    call la_dlasr('R','V','F',p,imax - imin + 1,work(iu1cs + imin - 1),work( &
                              iu1sn + imin - 1),u1(1,imin),ldu1)
                 else
                    call la_dlasr('L','V','F',imax - imin + 1,p,work(iu1cs + imin - 1),work( &
                              iu1sn + imin - 1),u1(imin,1),ldu1)
                 end if
              end if
              if (wantu2) then
                 if (colmajor) then
                    call la_dlasr('R','V','F',m - p,imax - imin + 1,work(iu2cs + imin - 1),work( &
                              iu2sn + imin - 1),u2(1,imin),ldu2)
                 else
                    call la_dlasr('L','V','F',imax - imin + 1,m - p,work(iu2cs + imin - 1),work( &
                              iu2sn + imin - 1),u2(imin,1),ldu2)
                 end if
              end if
              if (wantv1t) then
                 if (colmajor) then
                    call la_dlasr('L','V','F',imax - imin + 1,q,work(iv1tcs + imin - 1),work( &
                              iv1tsn + imin - 1),v1t(imin,1),ldv1t)
                 else
                    call la_dlasr('R','V','F',q,imax - imin + 1,work(iv1tcs + imin - 1),work( &
                              iv1tsn + imin - 1),v1t(1,imin),ldv1t)
                 end if
              end if
              if (wantv2t) then
                 if (colmajor) then
                    call la_dlasr('L','V','F',imax - imin + 1,m - q,work(iv2tcs + imin - 1),work( &
                              iv2tsn + imin - 1),v2t(imin,1),ldv2t)
                 else
                    call la_dlasr('R','V','F',m - q,imax - imin + 1,work(iv2tcs + imin - 1),work( &
                              iv2tsn + imin - 1),v2t(1,imin),ldv2t)
                 end if
              end if
              ! fix signs on b11(imax-1,imax) and b21(imax-1,imax)
              if (b11e(imax - 1) + b21e(imax - 1) > 0) then
                 b11d(imax) = -b11d(imax)
                 b21d(imax) = -b21d(imax)
                 if (wantv1t) then
                    if (colmajor) then
                       call la_dscal(q,negone,v1t(imax,1),ldv1t)
                    else
                       call la_dscal(q,negone,v1t(1,imax),1)
                    end if
                 end if
              end if
              ! compute theta(imax)
              x1 = cos(phi(imax - 1))*b11d(imax) + sin(phi(imax - 1))*b12e(imax - 1)
              y1 = cos(phi(imax - 1))*b21d(imax) + sin(phi(imax - 1))*b22e(imax - 1)
              theta(imax) = atan2(abs(y1),abs(x1))
              ! fix signs on b11(imax,imax), b12(imax,imax-1), b21(imax,imax),
              ! and b22(imax,imax-1)
              if (b11d(imax) + b12e(imax - 1) < 0) then
                 b12d(imax) = -b12d(imax)
                 if (wantu1) then
                    if (colmajor) then
                       call la_dscal(p,negone,u1(1,imax),1)
                    else
                       call la_dscal(p,negone,u1(imax,1),ldu1)
                    end if
                 end if
              end if
              if (b21d(imax) + b22e(imax - 1) > 0) then
                 b22d(imax) = -b22d(imax)
                 if (wantu2) then
                    if (colmajor) then
                       call la_dscal(m - p,negone,u2(1,imax),1)
                    else
                       call la_dscal(m - p,negone,u2(imax,1),ldu2)
                    end if
                 end if
              end if
              ! fix signs on b12(imax,imax) and b22(imax,imax)
              if (b12d(imax) + b22d(imax) < 0) then
                 if (wantv2t) then
                    if (colmajor) then
                       call la_dscal(m - q,negone,v2t(imax,1),ldv2t)
                    else
                       call la_dscal(m - q,negone,v2t(1,imax),1)
                    end if
                 end if
              end if
              ! test for negligible sines or cosines
              do i = imin,imax
                 if (theta(i) < thresh) then
                    theta(i) = zero
                 else if (theta(i) > piover2 - thresh) then
                    theta(i) = piover2
                 end if
              end do
              do i = imin,imax - 1
                 if (phi(i) < thresh) then
                    phi(i) = zero
                 else if (phi(i) > piover2 - thresh) then
                    phi(i) = piover2
                 end if
              end do
              ! deflate
              if (imax > 1) then
                 do while (phi(imax - 1) == zero)
                    imax = imax - 1
                    if (imax <= 1) exit
                 end do
              end if
              if (imin > imax - 1) imin = imax - 1
              if (imin > 1) then
                 do while (phi(imin - 1) /= zero)
                     imin = imin - 1
                     if (imin <= 1) exit
                 end do
              end if
              ! repeat main iteration loop
           end do
           ! postprocessing: order theta from least to greatest
           do i = 1,q
              mini = i
              thetamin = theta(i)
              do j = i + 1,q
                 if (theta(j) < thetamin) then
                    mini = j
                    thetamin = theta(j)
                 end if
              end do
              if (mini /= i) then
                 theta(mini) = theta(i)
                 theta(i) = thetamin
                 if (colmajor) then
                    if (wantu1) call la_dswap(p,u1(1,i),1,u1(1,mini),1)
                    if (wantu2) call la_dswap(m - p,u2(1,i),1,u2(1,mini),1)
                    if (wantv1t) call la_dswap(q,v1t(i,1),ldv1t,v1t(mini,1),ldv1t)

                    if (wantv2t) call la_dswap(m - q,v2t(i,1),ldv2t,v2t(mini,1),ldv2t)

                 else
                    if (wantu1) call la_dswap(p,u1(i,1),ldu1,u1(mini,1),ldu1)
                    if (wantu2) call la_dswap(m - p,u2(i,1),ldu2,u2(mini,1),ldu2)
                    if (wantv1t) call la_dswap(q,v1t(1,i),1,v1t(1,mini),1)
                    if (wantv2t) call la_dswap(m - q,v2t(1,i),1,v2t(1,mini),1)
                 end if
              end if
           end do
           return
     end subroutine la_dbbcsd
     !> QBBCSD: computes the CS decomposition of an orthogonal matrix in
     !> bidiagonal-block form,
     !> [ B11 | B12 0  0 ]
     !> [  0  |  0 -I  0 ]
     !> X = [----------------]
     !> [ B21 | B22 0  0 ]
     !> [  0  |  0  0  I ]
     !> [  C | -S  0  0 ]
     !> [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T
     !> = [---------] [---------------] [---------]   .
     !> [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
     !> [  0 |  0  0  I ]
     !> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
     !> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
     !> transposed and/or permuted. This can be done in constant time using
     !> the TRANS and SIGNS options. See QORCSD for details.)
     !> The bidiagonal matrices B11, B12, B21, and B22 are represented
     !> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
     !> The orthogonal matrices U1, U2, V1T, and V2T are input/output.
     !> The input matrices are pre- or post-multiplied by the appropriate
     !> singular vector matrices.

     pure subroutine la_qbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,phi,u1, &
     ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,b11d,b11e,b12d,b12e,b21d,b21e,b22d,b22e,work, &
               lwork,info)
        use la_constants_qp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,lwork,m,p,q
           ! Array Arguments
           real(qp),intent(out) :: b11d(*),b11e(*),b12d(*),b12e(*),b21d(*),b21e(*),b22d(*), &
                      b22e(*),work(*)
           real(qp),intent(inout) :: phi(*),theta(*)
           real(qp),intent(inout) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*)
        ! ===================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 6
           real(qp),parameter :: hundred = 100.0_qp
           real(qp),parameter :: meighth = -0.125_qp
           real(qp),parameter :: piover2 = 1.57079632679489661923132169163975144210_qp

           ! Local Scalars
           logical(lk) :: colmajor,lquery,restart11,restart12,restart21,restart22,wantu1, &
                     wantu2,wantv1t,wantv2t
           integer(ilp) :: i,imin,imax,iter,iu1cs,iu1sn,iu2cs,iu2sn,iv1tcs,iv1tsn, &
                     iv2tcs,iv2tsn,j,lworkmin,lworkopt,maxit,mini
           real(qp) :: b11bulge,b12bulge,b21bulge,b22bulge,dummy,eps,mu,nu,r,sigma11, &
                     sigma21,temp,thetamax,thetamin,thresh,tol,tolmul,unfl,x1,x2,y1,y2
           ! Intrinsic Functions
           intrinsic :: abs,atan2,cos,max,min,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           if (m < 0) then
              info = -6
           else if (p < 0 .or. p > m) then
              info = -7
           else if (q < 0 .or. q > m) then
              info = -8
           else if (q > p .or. q > m - p .or. q > m - q) then
              info = -8
           else if (wantu1 .and. ldu1 < p) then
              info = -12
           else if (wantu2 .and. ldu2 < m - p) then
              info = -14
           else if (wantv1t .and. ldv1t < q) then
              info = -16
           else if (wantv2t .and. ldv2t < m - q) then
              info = -18
           end if
           ! quick return if q = 0
           if (info == 0 .and. q == 0) then
              lworkmin = 1
              work(1) = lworkmin
              return
           end if
           ! compute workspace
           if (info == 0) then
              iu1cs = 1
              iu1sn = iu1cs + q
              iu2cs = iu1sn + q
              iu2sn = iu2cs + q
              iv1tcs = iu2sn + q
              iv1tsn = iv1tcs + q
              iv2tcs = iv1tsn + q
              iv2tsn = iv2tcs + q
              lworkopt = iv2tsn + q - 1
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -28
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QBBCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! get machine constants
           eps = la_qlamch('EPSILON')
           unfl = la_qlamch('SAFE MINIMUM')
           tolmul = max(ten,min(hundred,eps**meighth))
           tol = tolmul*eps
           thresh = max(tol,maxitr*q*q*unfl)
           ! test for negligible sines or cosines
           do i = 1,q
              if (theta(i) < thresh) then
                 theta(i) = zero
              else if (theta(i) > piover2 - thresh) then
                 theta(i) = piover2
              end if
           end do
           do i = 1,q - 1
              if (phi(i) < thresh) then
                 phi(i) = zero
              else if (phi(i) > piover2 - thresh) then
                 phi(i) = piover2
              end if
           end do
           ! initial deflation
           imax = q
           do while (imax > 1)
              if (phi(imax - 1) /= zero) then
                 exit
              end if
              imax = imax - 1
           end do
           imin = imax - 1
           if (imin > 1) then
              do while (phi(imin - 1) /= zero)
                 imin = imin - 1
                 if (imin <= 1) exit
              end do
           end if
           ! initialize iteration counter
           maxit = maxitr*q*q
           iter = 0
           ! begin main iteration loop
           do while (imax > 1)
              ! compute the matrix entries
              b11d(imin) = cos(theta(imin))
              b21d(imin) = -sin(theta(imin))
              do i = imin,imax - 1
                 b11e(i) = -sin(theta(i))*sin(phi(i))
                 b11d(i + 1) = cos(theta(i + 1))*cos(phi(i))
                 b12d(i) = sin(theta(i))*cos(phi(i))
                 b12e(i) = cos(theta(i + 1))*sin(phi(i))
                 b21e(i) = -cos(theta(i))*sin(phi(i))
                 b21d(i + 1) = -sin(theta(i + 1))*cos(phi(i))
                 b22d(i) = cos(theta(i))*cos(phi(i))
                 b22e(i) = -sin(theta(i + 1))*sin(phi(i))
              end do
              b12d(imax) = sin(theta(imax))
              b22d(imax) = cos(theta(imax))
              ! abort if not converging; otherwise, increment iter
              if (iter > maxit) then
                 info = 0
                 do i = 1,q
                    if (phi(i) /= zero) info = info + 1
                 end do
                 return
              end if
              iter = iter + imax - imin
              ! compute shifts
              thetamax = theta(imin)
              thetamin = theta(imin)
              do i = imin + 1,imax
                 if (theta(i) > thetamax) thetamax = theta(i)
                 if (theta(i) < thetamin) thetamin = theta(i)
              end do
              if (thetamax > piover2 - thresh) then
                 ! zero on diagonals of b11 and b22; induce deflation with a
                 ! zero shift
                 mu = zero
                 nu = one
              else if (thetamin < thresh) then
                 ! zero on diagonals of b12 and b22; induce deflation with a
                 ! zero shift
                 mu = one
                 nu = zero
              else
                 ! compute shifts for b11 and b21 and use the lesser
                 call la_qlas2(b11d(imax - 1),b11e(imax - 1),b11d(imax),sigma11,dummy)

                 call la_qlas2(b21d(imax - 1),b21e(imax - 1),b21d(imax),sigma21,dummy)

                 if (sigma11 <= sigma21) then
                    mu = sigma11
                    nu = sqrt(one - mu**2)
                    if (mu < thresh) then
                       mu = zero
                       nu = one
                    end if
                 else
                    nu = sigma21
                    mu = sqrt(1.0_qp - nu**2)
                    if (nu < thresh) then
                       mu = one
                       nu = zero
                    end if
                 end if
              end if
              ! rotate to produce bulges in b11 and b21
              if (mu <= nu) then
                 call la_qlartgs(b11d(imin),b11e(imin),mu,work(iv1tcs + imin - 1),work(iv1tsn + &
                           imin - 1))
              else
                 call la_qlartgs(b21d(imin),b21e(imin),nu,work(iv1tcs + imin - 1),work(iv1tsn + &
                           imin - 1))
              end if
              temp = work(iv1tcs + imin - 1)*b11d(imin) + work(iv1tsn + imin - 1)*b11e(imin)
              b11e(imin) = work(iv1tcs + imin - 1)*b11e(imin) - work(iv1tsn + imin - 1)*b11d(imin)
              b11d(imin) = temp
              b11bulge = work(iv1tsn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = work(iv1tcs + imin - 1)*b11d(imin + 1)
              temp = work(iv1tcs + imin - 1)*b21d(imin) + work(iv1tsn + imin - 1)*b21e(imin)
              b21e(imin) = work(iv1tcs + imin - 1)*b21e(imin) - work(iv1tsn + imin - 1)*b21d(imin)
              b21d(imin) = temp
              b21bulge = work(iv1tsn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = work(iv1tcs + imin - 1)*b21d(imin + 1)
              ! compute theta(imin)
              theta(imin) = atan2(sqrt(b21d(imin)**2 + b21bulge**2),sqrt(b11d(imin)**2 + &
                        b11bulge**2))
              ! chase the bulges in b11(imin+1,imin) and b21(imin+1,imin)
              if (b11d(imin)**2 + b11bulge**2 > thresh**2) then
                 call la_qlartgp(b11bulge,b11d(imin),work(iu1sn + imin - 1),work(iu1cs + imin - 1), &
                            r)
              else if (mu <= nu) then
                 call la_qlartgs(b11e(imin),b11d(imin + 1),mu,work(iu1cs + imin - 1),work( &
                           iu1sn + imin - 1))
              else
                 call la_qlartgs(b12d(imin),b12e(imin),nu,work(iu1cs + imin - 1),work( &
                           iu1sn + imin - 1))
              end if
              if (b21d(imin)**2 + b21bulge**2 > thresh**2) then
                 call la_qlartgp(b21bulge,b21d(imin),work(iu2sn + imin - 1),work(iu2cs + imin - 1), &
                            r)
              else if (nu < mu) then
                 call la_qlartgs(b21e(imin),b21d(imin + 1),nu,work(iu2cs + imin - 1),work( &
                           iu2sn + imin - 1))
              else
                 call la_qlartgs(b22d(imin),b22e(imin),mu,work(iu2cs + imin - 1),work(iu2sn + &
                           imin - 1))
              end if
              work(iu2cs + imin - 1) = -work(iu2cs + imin - 1)
              work(iu2sn + imin - 1) = -work(iu2sn + imin - 1)
              temp = work(iu1cs + imin - 1)*b11e(imin) + work(iu1sn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = work(iu1cs + imin - 1)*b11d(imin + 1) - work(iu1sn + imin - 1)*b11e(imin)

              b11e(imin) = temp
              if (imax > imin + 1) then
                 b11bulge = work(iu1sn + imin - 1)*b11e(imin + 1)
                 b11e(imin + 1) = work(iu1cs + imin - 1)*b11e(imin + 1)
              end if
              temp = work(iu1cs + imin - 1)*b12d(imin) + work(iu1sn + imin - 1)*b12e(imin)
              b12e(imin) = work(iu1cs + imin - 1)*b12e(imin) - work(iu1sn + imin - 1)*b12d(imin)
              b12d(imin) = temp
              b12bulge = work(iu1sn + imin - 1)*b12d(imin + 1)
              b12d(imin + 1) = work(iu1cs + imin - 1)*b12d(imin + 1)
              temp = work(iu2cs + imin - 1)*b21e(imin) + work(iu2sn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = work(iu2cs + imin - 1)*b21d(imin + 1) - work(iu2sn + imin - 1)*b21e(imin)

              b21e(imin) = temp
              if (imax > imin + 1) then
                 b21bulge = work(iu2sn + imin - 1)*b21e(imin + 1)
                 b21e(imin + 1) = work(iu2cs + imin - 1)*b21e(imin + 1)
              end if
              temp = work(iu2cs + imin - 1)*b22d(imin) + work(iu2sn + imin - 1)*b22e(imin)
              b22e(imin) = work(iu2cs + imin - 1)*b22e(imin) - work(iu2sn + imin - 1)*b22d(imin)
              b22d(imin) = temp
              b22bulge = work(iu2sn + imin - 1)*b22d(imin + 1)
              b22d(imin + 1) = work(iu2cs + imin - 1)*b22d(imin + 1)
              ! inner loop: chase bulges from b11(imin,imin+2),
              ! b12(imin,imin+1), b21(imin,imin+2), and b22(imin,imin+1) to
              ! bottom-right
              do i = imin + 1,imax - 1
                 ! compute phi(i-1)
                 x1 = sin(theta(i - 1))*b11e(i - 1) + cos(theta(i - 1))*b21e(i - 1)
                 x2 = sin(theta(i - 1))*b11bulge + cos(theta(i - 1))*b21bulge
                 y1 = sin(theta(i - 1))*b12d(i - 1) + cos(theta(i - 1))*b22d(i - 1)
                 y2 = sin(theta(i - 1))*b12bulge + cos(theta(i - 1))*b22bulge
                 phi(i - 1) = atan2(sqrt(x1**2 + x2**2),sqrt(y1**2 + y2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11e(i - 1)**2 + b11bulge**2 <= thresh**2
                 restart21 = b21e(i - 1)**2 + b21bulge**2 <= thresh**2
                 restart12 = b12d(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart22 = b22d(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i-1,i+1), b12(i-1,i),
                 ! b21(i-1,i+1), and b22(i-1,i). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart21) then
                    call la_qlartgp(x2,x1,work(iv1tsn + i - 1),work(iv1tcs + i - 1),r)
                 else if (.not. restart11 .and. restart21) then
                    call la_qlartgp(b11bulge,b11e(i - 1),work(iv1tsn + i - 1),work(iv1tcs + i - 1), &
                              r)
                 else if (restart11 .and. .not. restart21) then
                    call la_qlartgp(b21bulge,b21e(i - 1),work(iv1tsn + i - 1),work(iv1tcs + i - 1), &
                              r)
                 else if (mu <= nu) then
                    call la_qlartgs(b11d(i),b11e(i),mu,work(iv1tcs + i - 1),work(iv1tsn + i - 1))

                 else
                    call la_qlartgs(b21d(i),b21e(i),nu,work(iv1tcs + i - 1),work(iv1tsn + i - 1))

                 end if
                 work(iv1tcs + i - 1) = -work(iv1tcs + i - 1)
                 work(iv1tsn + i - 1) = -work(iv1tsn + i - 1)
                 if (.not. restart12 .and. .not. restart22) then
                    call la_qlartgp(y2,y1,work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - 1),r)

                 else if (.not. restart12 .and. restart22) then
                    call la_qlartgp(b12bulge,b12d(i - 1),work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - &
                              1),r)
                 else if (restart12 .and. .not. restart22) then
                    call la_qlartgp(b22bulge,b22d(i - 1),work(iv2tsn + i - 1 - 1),work(iv2tcs + i - 1 - &
                              1),r)
                 else if (nu < mu) then
                    call la_qlartgs(b12e(i - 1),b12d(i),nu,work(iv2tcs + i - 1 - 1),work(iv2tsn + i - &
                              1 - 1))
                 else
                    call la_qlartgs(b22e(i - 1),b22d(i),mu,work(iv2tcs + i - 1 - 1),work(iv2tsn + i - &
                              1 - 1))
                 end if
                 temp = work(iv1tcs + i - 1)*b11d(i) + work(iv1tsn + i - 1)*b11e(i)
                 b11e(i) = work(iv1tcs + i - 1)*b11e(i) - work(iv1tsn + i - 1)*b11d(i)
                 b11d(i) = temp
                 b11bulge = work(iv1tsn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = work(iv1tcs + i - 1)*b11d(i + 1)
                 temp = work(iv1tcs + i - 1)*b21d(i) + work(iv1tsn + i - 1)*b21e(i)
                 b21e(i) = work(iv1tcs + i - 1)*b21e(i) - work(iv1tsn + i - 1)*b21d(i)
                 b21d(i) = temp
                 b21bulge = work(iv1tsn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = work(iv1tcs + i - 1)*b21d(i + 1)
                 temp = work(iv2tcs + i - 1 - 1)*b12e(i - 1) + work(iv2tsn + i - 1 - 1)*b12d(i)
                 b12d(i) = work(iv2tcs + i - 1 - 1)*b12d(i) - work(iv2tsn + i - 1 - 1)*b12e(i - 1)
                 b12e(i - 1) = temp
                 b12bulge = work(iv2tsn + i - 1 - 1)*b12e(i)
                 b12e(i) = work(iv2tcs + i - 1 - 1)*b12e(i)
                 temp = work(iv2tcs + i - 1 - 1)*b22e(i - 1) + work(iv2tsn + i - 1 - 1)*b22d(i)
                 b22d(i) = work(iv2tcs + i - 1 - 1)*b22d(i) - work(iv2tsn + i - 1 - 1)*b22e(i - 1)
                 b22e(i - 1) = temp
                 b22bulge = work(iv2tsn + i - 1 - 1)*b22e(i)
                 b22e(i) = work(iv2tcs + i - 1 - 1)*b22e(i)
                 ! compute theta(i)
                 x1 = cos(phi(i - 1))*b11d(i) + sin(phi(i - 1))*b12e(i - 1)
                 x2 = cos(phi(i - 1))*b11bulge + sin(phi(i - 1))*b12bulge
                 y1 = cos(phi(i - 1))*b21d(i) + sin(phi(i - 1))*b22e(i - 1)
                 y2 = cos(phi(i - 1))*b21bulge + sin(phi(i - 1))*b22bulge
                 theta(i) = atan2(sqrt(y1**2 + y2**2),sqrt(x1**2 + x2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11d(i)**2 + b11bulge**2 <= thresh**2
                 restart12 = b12e(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart21 = b21d(i)**2 + b21bulge**2 <= thresh**2
                 restart22 = b22e(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i+1,i), b12(i+1,i-1),
                 ! b21(i+1,i), and b22(i+1,i-1). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart12) then
                    call la_qlartgp(x2,x1,work(iu1sn + i - 1),work(iu1cs + i - 1),r)
                 else if (.not. restart11 .and. restart12) then
                    call la_qlartgp(b11bulge,b11d(i),work(iu1sn + i - 1),work(iu1cs + i - 1),r)

                 else if (restart11 .and. .not. restart12) then
                    call la_qlartgp(b12bulge,b12e(i - 1),work(iu1sn + i - 1),work(iu1cs + i - 1),r)

                 else if (mu <= nu) then
                    call la_qlartgs(b11e(i),b11d(i + 1),mu,work(iu1cs + i - 1),work(iu1sn + i - 1))

                 else
                    call la_qlartgs(b12d(i),b12e(i),nu,work(iu1cs + i - 1),work(iu1sn + i - 1))

                 end if
                 if (.not. restart21 .and. .not. restart22) then
                    call la_qlartgp(y2,y1,work(iu2sn + i - 1),work(iu2cs + i - 1),r)
                 else if (.not. restart21 .and. restart22) then
                    call la_qlartgp(b21bulge,b21d(i),work(iu2sn + i - 1),work(iu2cs + i - 1),r)

                 else if (restart21 .and. .not. restart22) then
                    call la_qlartgp(b22bulge,b22e(i - 1),work(iu2sn + i - 1),work(iu2cs + i - 1),r)

                 else if (nu < mu) then
                    call la_qlartgs(b21e(i),b21e(i + 1),nu,work(iu2cs + i - 1),work(iu2sn + i - 1))

                 else
                    call la_qlartgs(b22d(i),b22e(i),mu,work(iu2cs + i - 1),work(iu2sn + i - 1))

                 end if
                 work(iu2cs + i - 1) = -work(iu2cs + i - 1)
                 work(iu2sn + i - 1) = -work(iu2sn + i - 1)
                 temp = work(iu1cs + i - 1)*b11e(i) + work(iu1sn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = work(iu1cs + i - 1)*b11d(i + 1) - work(iu1sn + i - 1)*b11e(i)
                 b11e(i) = temp
                 if (i < imax - 1) then
                    b11bulge = work(iu1sn + i - 1)*b11e(i + 1)
                    b11e(i + 1) = work(iu1cs + i - 1)*b11e(i + 1)
                 end if
                 temp = work(iu2cs + i - 1)*b21e(i) + work(iu2sn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = work(iu2cs + i - 1)*b21d(i + 1) - work(iu2sn + i - 1)*b21e(i)
                 b21e(i) = temp
                 if (i < imax - 1) then
                    b21bulge = work(iu2sn + i - 1)*b21e(i + 1)
                    b21e(i + 1) = work(iu2cs + i - 1)*b21e(i + 1)
                 end if
                 temp = work(iu1cs + i - 1)*b12d(i) + work(iu1sn + i - 1)*b12e(i)
                 b12e(i) = work(iu1cs + i - 1)*b12e(i) - work(iu1sn + i - 1)*b12d(i)
                 b12d(i) = temp
                 b12bulge = work(iu1sn + i - 1)*b12d(i + 1)
                 b12d(i + 1) = work(iu1cs + i - 1)*b12d(i + 1)
                 temp = work(iu2cs + i - 1)*b22d(i) + work(iu2sn + i - 1)*b22e(i)
                 b22e(i) = work(iu2cs + i - 1)*b22e(i) - work(iu2sn + i - 1)*b22d(i)
                 b22d(i) = temp
                 b22bulge = work(iu2sn + i - 1)*b22d(i + 1)
                 b22d(i + 1) = work(iu2cs + i - 1)*b22d(i + 1)
              end do
              ! compute phi(imax-1)
              x1 = sin(theta(imax - 1))*b11e(imax - 1) + cos(theta(imax - 1))*b21e(imax - 1)
              y1 = sin(theta(imax - 1))*b12d(imax - 1) + cos(theta(imax - 1))*b22d(imax - 1)
              y2 = sin(theta(imax - 1))*b12bulge + cos(theta(imax - 1))*b22bulge
              phi(imax - 1) = atan2(abs(x1),sqrt(y1**2 + y2**2))
              ! chase bulges from b12(imax-1,imax) and b22(imax-1,imax)
              restart12 = b12d(imax - 1)**2 + b12bulge**2 <= thresh**2
              restart22 = b22d(imax - 1)**2 + b22bulge**2 <= thresh**2
              if (.not. restart12 .and. .not. restart22) then
                 call la_qlartgp(y2,y1,work(iv2tsn + imax - 1 - 1),work(iv2tcs + imax - 1 - 1),r)

              else if (.not. restart12 .and. restart22) then
                 call la_qlartgp(b12bulge,b12d(imax - 1),work(iv2tsn + imax - 1 - 1),work(iv2tcs + &
                           imax - 1 - 1),r)
              else if (restart12 .and. .not. restart22) then
                 call la_qlartgp(b22bulge,b22d(imax - 1),work(iv2tsn + imax - 1 - 1),work(iv2tcs + &
                           imax - 1 - 1),r)
              else if (nu < mu) then
                 call la_qlartgs(b12e(imax - 1),b12d(imax),nu,work(iv2tcs + imax - 1 - 1),work( &
                           iv2tsn + imax - 1 - 1))
              else
                 call la_qlartgs(b22e(imax - 1),b22d(imax),mu,work(iv2tcs + imax - 1 - 1),work( &
                           iv2tsn + imax - 1 - 1))
              end if
              temp = work(iv2tcs + imax - 1 - 1)*b12e(imax - 1) + work(iv2tsn + imax - 1 - 1)*b12d(imax)
              b12d(imax) = work(iv2tcs + imax - 1 - 1)*b12d(imax) - work(iv2tsn + imax - 1 - 1)*b12e(imax - 1)

              b12e(imax - 1) = temp
              temp = work(iv2tcs + imax - 1 - 1)*b22e(imax - 1) + work(iv2tsn + imax - 1 - 1)*b22d(imax)
              b22d(imax) = work(iv2tcs + imax - 1 - 1)*b22d(imax) - work(iv2tsn + imax - 1 - 1)*b22e(imax - 1)

              b22e(imax - 1) = temp
              ! update singular vectors
              if (wantu1) then
                 if (colmajor) then
                    call la_qlasr('R','V','F',p,imax - imin + 1,work(iu1cs + imin - 1),work( &
                              iu1sn + imin - 1),u1(1,imin),ldu1)
                 else
                    call la_qlasr('L','V','F',imax - imin + 1,p,work(iu1cs + imin - 1),work( &
                              iu1sn + imin - 1),u1(imin,1),ldu1)
                 end if
              end if
              if (wantu2) then
                 if (colmajor) then
                    call la_qlasr('R','V','F',m - p,imax - imin + 1,work(iu2cs + imin - 1),work( &
                              iu2sn + imin - 1),u2(1,imin),ldu2)
                 else
                    call la_qlasr('L','V','F',imax - imin + 1,m - p,work(iu2cs + imin - 1),work( &
                              iu2sn + imin - 1),u2(imin,1),ldu2)
                 end if
              end if
              if (wantv1t) then
                 if (colmajor) then
                    call la_qlasr('L','V','F',imax - imin + 1,q,work(iv1tcs + imin - 1),work( &
                              iv1tsn + imin - 1),v1t(imin,1),ldv1t)
                 else
                    call la_qlasr('R','V','F',q,imax - imin + 1,work(iv1tcs + imin - 1),work( &
                              iv1tsn + imin - 1),v1t(1,imin),ldv1t)
                 end if
              end if
              if (wantv2t) then
                 if (colmajor) then
                    call la_qlasr('L','V','F',imax - imin + 1,m - q,work(iv2tcs + imin - 1),work( &
                              iv2tsn + imin - 1),v2t(imin,1),ldv2t)
                 else
                    call la_qlasr('R','V','F',m - q,imax - imin + 1,work(iv2tcs + imin - 1),work( &
                              iv2tsn + imin - 1),v2t(1,imin),ldv2t)
                 end if
              end if
              ! fix signs on b11(imax-1,imax) and b21(imax-1,imax)
              if (b11e(imax - 1) + b21e(imax - 1) > 0) then
                 b11d(imax) = -b11d(imax)
                 b21d(imax) = -b21d(imax)
                 if (wantv1t) then
                    if (colmajor) then
                       call la_qscal(q,negone,v1t(imax,1),ldv1t)
                    else
                       call la_qscal(q,negone,v1t(1,imax),1)
                    end if
                 end if
              end if
              ! compute theta(imax)
              x1 = cos(phi(imax - 1))*b11d(imax) + sin(phi(imax - 1))*b12e(imax - 1)
              y1 = cos(phi(imax - 1))*b21d(imax) + sin(phi(imax - 1))*b22e(imax - 1)
              theta(imax) = atan2(abs(y1),abs(x1))
              ! fix signs on b11(imax,imax), b12(imax,imax-1), b21(imax,imax),
              ! and b22(imax,imax-1)
              if (b11d(imax) + b12e(imax - 1) < 0) then
                 b12d(imax) = -b12d(imax)
                 if (wantu1) then
                    if (colmajor) then
                       call la_qscal(p,negone,u1(1,imax),1)
                    else
                       call la_qscal(p,negone,u1(imax,1),ldu1)
                    end if
                 end if
              end if
              if (b21d(imax) + b22e(imax - 1) > 0) then
                 b22d(imax) = -b22d(imax)
                 if (wantu2) then
                    if (colmajor) then
                       call la_qscal(m - p,negone,u2(1,imax),1)
                    else
                       call la_qscal(m - p,negone,u2(imax,1),ldu2)
                    end if
                 end if
              end if
              ! fix signs on b12(imax,imax) and b22(imax,imax)
              if (b12d(imax) + b22d(imax) < 0) then
                 if (wantv2t) then
                    if (colmajor) then
                       call la_qscal(m - q,negone,v2t(imax,1),ldv2t)
                    else
                       call la_qscal(m - q,negone,v2t(1,imax),1)
                    end if
                 end if
              end if
              ! test for negligible sines or cosines
              do i = imin,imax
                 if (theta(i) < thresh) then
                    theta(i) = zero
                 else if (theta(i) > piover2 - thresh) then
                    theta(i) = piover2
                 end if
              end do
              do i = imin,imax - 1
                 if (phi(i) < thresh) then
                    phi(i) = zero
                 else if (phi(i) > piover2 - thresh) then
                    phi(i) = piover2
                 end if
              end do
              ! deflate
              if (imax > 1) then
                 do while (phi(imax - 1) == zero)
                    imax = imax - 1
                    if (imax <= 1) exit
                 end do
              end if
              if (imin > imax - 1) imin = imax - 1
              if (imin > 1) then
                 do while (phi(imin - 1) /= zero)
                     imin = imin - 1
                     if (imin <= 1) exit
                 end do
              end if
              ! repeat main iteration loop
           end do
           ! postprocessing: order theta from least to greatest
           do i = 1,q
              mini = i
              thetamin = theta(i)
              do j = i + 1,q
                 if (theta(j) < thetamin) then
                    mini = j
                    thetamin = theta(j)
                 end if
              end do
              if (mini /= i) then
                 theta(mini) = theta(i)
                 theta(i) = thetamin
                 if (colmajor) then
                    if (wantu1) call la_qswap(p,u1(1,i),1,u1(1,mini),1)
                    if (wantu2) call la_qswap(m - p,u2(1,i),1,u2(1,mini),1)
                    if (wantv1t) call la_qswap(q,v1t(i,1),ldv1t,v1t(mini,1),ldv1t)

                    if (wantv2t) call la_qswap(m - q,v2t(i,1),ldv2t,v2t(mini,1),ldv2t)

                 else
                    if (wantu1) call la_qswap(p,u1(i,1),ldu1,u1(mini,1),ldu1)
                    if (wantu2) call la_qswap(m - p,u2(i,1),ldu2,u2(mini,1),ldu2)
                    if (wantv1t) call la_qswap(q,v1t(1,i),1,v1t(1,mini),1)
                    if (wantv2t) call la_qswap(m - q,v2t(1,i),1,v2t(1,mini),1)
                 end if
              end if
           end do
           return
     end subroutine la_qbbcsd

     !> SORBDB: simultaneously bidiagonalizes the blocks of an M-by-M
     !> partitioned orthogonal matrix X:
     !> [ B11 | B12 0  0 ]
     !> [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T
     !> X = [-----------] = [---------] [----------------] [---------]   .
     !> [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
     !> [  0  |  0  0  I ]
     !> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     !> not the case, then X must be transposed and/or permuted. This can be
     !> done in constant time using the TRANS and SIGNS options. See SORCSD
     !> for details.)
     !> The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     !> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     !> represented implicitly by Householder vectors.
     !> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_sorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
               ldx22,theta,phi,taup1,taup2,tauq1,tauq2,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldx11,ldx12,ldx21,ldx22,lwork,m,p,q
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           real(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),tauq2(*),work(*)
           real(sp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ====================================================================
           ! Parameters
           real(sp),parameter :: realone = 1.0_sp

           ! Local Scalars
           logical(lk) :: colmajor,lquery
           integer(ilp) :: i,lworkmin,lworkopt
           real(sp) :: z1,z2,z3,z4
           ! Intrinsic Functions
           intrinsic :: atan2,cos,max,sin
           ! Executable Statements
           ! test input arguments
           info = 0
           colmajor = .not. la_lsame(trans,'T')
           if (.not. la_lsame(signs,'O')) then
              z1 = realone
              z2 = realone
              z3 = realone
              z4 = realone
           else
              z1 = realone
              z2 = -realone
              z3 = realone
              z4 = -realone
           end if
           lquery = lwork == -1
           if (m < 0) then
              info = -3
           else if (p < 0 .or. p > m) then
              info = -4
           else if (q < 0 .or. q > p .or. q > m - p .or. q > m - q) then
              info = -5
           else if (colmajor .and. ldx11 < max(1,p)) then
              info = -7
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
              info = -7
           else if (colmajor .and. ldx12 < max(1,p)) then
              info = -9
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
              info = -9
           else if (colmajor .and. ldx21 < max(1,m - p)) then
              info = -11
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
              info = -11
           else if (colmajor .and. ldx22 < max(1,m - p)) then
              info = -13
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
              info = -13
           end if
           ! compute workspace
           if (info == 0) then
              lworkopt = m - q
              lworkmin = m - q
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORBDB',-info)
              return
           else if (lquery) then
              return
           end if
           ! handle column-major and row-major separately
           if (colmajor) then
              ! reduce columns 1, ..., q of x11, x12, x21, and x22
              do i = 1,q
                 if (i == 1) then
                    call la_sscal(p - i + 1,z1,x11(i,i),1)
                 else
                    call la_sscal(p - i + 1,z1*cos(phi(i - 1)),x11(i,i),1)
                    call la_saxpy(p - i + 1,-z1*z3*z4*sin(phi(i - 1)),x12(i,i - 1),1,x11(i,i),1)

                 end if
                 if (i == 1) then
                    call la_sscal(m - p - i + 1,z2,x21(i,i),1)
                 else
                    call la_sscal(m - p - i + 1,z2*cos(phi(i - 1)),x21(i,i),1)
                    call la_saxpy(m - p - i + 1,-z2*z3*z4*sin(phi(i - 1)),x22(i,i - 1),1,x21(i,i), &
                              1)
                 end if
                 theta(i) = atan2(la_snrm2(m - p - i + 1,x21(i,i),1),la_snrm2(p - i + 1,x11( &
                           i,i),1))
                 if (p > i) then
                    call la_slarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
                 else if (p == i) then
                    call la_slarfgp(p - i + 1,x11(i,i),x11(i,i),1,taup1(i))
                 end if
                 x11(i,i) = one
                 if (m - p > i) then
                    call la_slarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
                 else if (m - p == i) then
                    call la_slarfgp(m - p - i + 1,x21(i,i),x21(i,i),1,taup2(i))
                 end if
                 x21(i,i) = one
                 if (q > i) then
                    call la_slarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11, &
                              work)
                 end if
                 if (m - q + 1 > i) then
                    call la_slarf('L',p - i + 1,m - q - i + 1,x11(i,i),1,taup1(i),x12(i,i),ldx12, &
                               work)
                 end if
                 if (q > i) then
                    call la_slarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21, &
                               work)
                 end if
                 if (m - q + 1 > i) then
                    call la_slarf('L',m - p - i + 1,m - q - i + 1,x21(i,i),1,taup2(i),x22(i,i), &
                              ldx22,work)
                 end if
                 if (i < q) then
                    call la_sscal(q - i,-z1*z3*sin(theta(i)),x11(i,i + 1),ldx11)
                    call la_saxpy(q - i,z2*z3*cos(theta(i)),x21(i,i + 1),ldx21,x11(i,i + 1), &
                              ldx11)
                 end if
                 call la_sscal(m - q - i + 1,-z1*z4*sin(theta(i)),x12(i,i),ldx12)
                 call la_saxpy(m - q - i + 1,z2*z4*cos(theta(i)),x22(i,i),ldx22,x12(i,i),ldx12 &
                           )
                 if (i < q) phi(i) = atan2(la_snrm2(q - i,x11(i,i + 1),ldx11),la_snrm2( &
                           m - q - i + 1,x12(i,i),ldx12))
                 if (i < q) then
                    if (q - i == 1) then
                       call la_slarfgp(q - i,x11(i,i + 1),x11(i,i + 1),ldx11,tauq1(i))
                    else
                       call la_slarfgp(q - i,x11(i,i + 1),x11(i,i + 2),ldx11,tauq1(i))
                    end if
                    x11(i,i + 1) = one
                 end if
                 if (q + i - 1 < m) then
                    if (m - q == i) then
                       call la_slarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                    else
                       call la_slarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))

                    end if
                 end if
                 x12(i,i) = one
                 if (i < q) then
                    call la_slarf('R',p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x11(i + 1,i + 1), &
                              ldx11,work)
                    call la_slarf('R',m - p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 if (p > i) then
                    call la_slarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p > i) then
                    call la_slarf('R',m - p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x22(i + 1,i), &
                              ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_sscal(m - q - i + 1,-z1*z4,x12(i,i),ldx12)
                 if (i >= m - q) then
                    call la_slarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                 else
                    call la_slarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))
                 end if
                 x12(i,i) = one
                 if (p > i) then
                    call la_slarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_slarf('R',m - p - q,m - q - i + 1,x12(i,i),ldx12,tauq2(i), &
                            x22(q + 1,i),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_sscal(m - p - q - i + 1,z2*z4,x22(q + i,p + i),ldx22)
                 if (i == m - p - q) then
                    call la_slarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i),ldx22,tauq2(p + i))

                 else
                    call la_slarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i + 1),ldx22,tauq2(p + i) &
                               )
                 end if
                 x22(q + i,p + i) = one
                 if (i < m - p - q) then
                    call la_slarf('R',m - p - q - i,m - p - q - i + 1,x22(q + i,p + i),ldx22,tauq2(p + i), &
                              x22(q + i + 1,p + i),ldx22,work)
                 end if
              end do
           else
              ! reduce columns 1, ..., q of x11, x12, x21, x22
              do i = 1,q
                 if (i == 1) then
                    call la_sscal(p - i + 1,z1,x11(i,i),ldx11)
                 else
                    call la_sscal(p - i + 1,z1*cos(phi(i - 1)),x11(i,i),ldx11)
                    call la_saxpy(p - i + 1,-z1*z3*z4*sin(phi(i - 1)),x12(i - 1,i),ldx12,x11(i,i), &
                               ldx11)
                 end if
                 if (i == 1) then
                    call la_sscal(m - p - i + 1,z2,x21(i,i),ldx21)
                 else
                    call la_sscal(m - p - i + 1,z2*cos(phi(i - 1)),x21(i,i),ldx21)
                    call la_saxpy(m - p - i + 1,-z2*z3*z4*sin(phi(i - 1)),x22(i - 1,i),ldx22,x21(i, &
                              i),ldx21)
                 end if
                 theta(i) = atan2(la_snrm2(m - p - i + 1,x21(i,i),ldx21),la_snrm2(p - i + 1, &
                           x11(i,i),ldx11))
                 call la_slarfgp(p - i + 1,x11(i,i),x11(i,i + 1),ldx11,taup1(i))
                 x11(i,i) = one
                 if (i == m - p) then
                    call la_slarfgp(m - p - i + 1,x21(i,i),x21(i,i),ldx21,taup2(i))
                 else
                    call la_slarfgp(m - p - i + 1,x21(i,i),x21(i,i + 1),ldx21,taup2(i))
                 end if
                 x21(i,i) = one
                 if (q > i) then
                    call la_slarf('R',q - i,p - i + 1,x11(i,i),ldx11,taup1(i),x11(i + 1,i), &
                              ldx11,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_slarf('R',m - q - i + 1,p - i + 1,x11(i,i),ldx11,taup1(i),x12(i,i), &
                              ldx12,work)
                 end if
                 if (q > i) then
                    call la_slarf('R',q - i,m - p - i + 1,x21(i,i),ldx21,taup2(i),x21(i + 1,i), &
                              ldx21,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_slarf('R',m - q - i + 1,m - p - i + 1,x21(i,i),ldx21,taup2(i),x22(i,i), &
                              ldx22,work)
                 end if
                 if (i < q) then
                    call la_sscal(q - i,-z1*z3*sin(theta(i)),x11(i + 1,i),1)
                    call la_saxpy(q - i,z2*z3*cos(theta(i)),x21(i + 1,i),1,x11(i + 1,i),1)

                 end if
                 call la_sscal(m - q - i + 1,-z1*z4*sin(theta(i)),x12(i,i),1)
                 call la_saxpy(m - q - i + 1,z2*z4*cos(theta(i)),x22(i,i),1,x12(i,i),1)

                 if (i < q) phi(i) = atan2(la_snrm2(q - i,x11(i + 1,i),1),la_snrm2(m - q - &
                           i + 1,x12(i,i),1))
                 if (i < q) then
                    if (q - i == 1) then
                       call la_slarfgp(q - i,x11(i + 1,i),x11(i + 1,i),1,tauq1(i))
                    else
                       call la_slarfgp(q - i,x11(i + 1,i),x11(i + 2,i),1,tauq1(i))
                    end if
                    x11(i + 1,i) = one
                 end if
                 if (m - q > i) then
                    call la_slarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 else
                    call la_slarfgp(m - q - i + 1,x12(i,i),x12(i,i),1,tauq2(i))
                 end if
                 x12(i,i) = one
                 if (i < q) then
                    call la_slarf('L',q - i,p - i,x11(i + 1,i),1,tauq1(i),x11(i + 1,i + 1),ldx11, &
                               work)
                    call la_slarf('L',q - i,m - p - i,x11(i + 1,i),1,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 call la_slarf('L',m - q - i + 1,p - i,x12(i,i),1,tauq2(i),x12(i,i + 1),ldx12, &
                           work)
                 if (m - p - i > 0) then
                    call la_slarf('L',m - q - i + 1,m - p - i,x12(i,i),1,tauq2(i),x22(i,i + 1), &
                              ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_sscal(m - q - i + 1,-z1*z4,x12(i,i),1)
                 call la_slarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = one
                 if (p > i) then
                    call la_slarf('L',m - q - i + 1,p - i,x12(i,i),1,tauq2(i),x12(i,i + 1),ldx12, &
                               work)
                 end if
                 if (m - p - q >= 1) call la_slarf('L',m - q - i + 1,m - p - q,x12(i,i),1,tauq2(i), &
                           x22(i,q + 1),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_sscal(m - p - q - i + 1,z2*z4,x22(p + i,q + i),1)
                 if (m - p - q == i) then
                    call la_slarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i,q + i),1,tauq2(p + i))

                    x22(p + i,q + i) = one
                 else
                    call la_slarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i + 1,q + i),1,tauq2(p + i))

                    x22(p + i,q + i) = one
                    call la_slarf('L',m - p - q - i + 1,m - p - q - i,x22(p + i,q + i),1,tauq2(p + i),x22(p + &
                              i,q + i + 1),ldx22,work)
                 end if
              end do
           end if
           return
     end subroutine la_sorbdb
     !> DORBDB: simultaneously bidiagonalizes the blocks of an M-by-M
     !> partitioned orthogonal matrix X:
     !> [ B11 | B12 0  0 ]
     !> [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T
     !> X = [-----------] = [---------] [----------------] [---------]   .
     !> [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
     !> [  0  |  0  0  I ]
     !> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     !> not the case, then X must be transposed and/or permuted. This can be
     !> done in constant time using the TRANS and SIGNS options. See DORCSD
     !> for details.)
     !> The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     !> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     !> represented implicitly by Householder vectors.
     !> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_dorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
               ldx22,theta,phi,taup1,taup2,tauq1,tauq2,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldx11,ldx12,ldx21,ldx22,lwork,m,p,q
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           real(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),tauq2(*),work(*)
           real(dp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ====================================================================
           ! Parameters
           real(dp),parameter :: realone = 1.0_dp

           ! Local Scalars
           logical(lk) :: colmajor,lquery
           integer(ilp) :: i,lworkmin,lworkopt
           real(dp) :: z1,z2,z3,z4
           ! Intrinsic Functions
           intrinsic :: atan2,cos,max,sin
           ! Executable Statements
           ! test input arguments
           info = 0
           colmajor = .not. la_lsame(trans,'T')
           if (.not. la_lsame(signs,'O')) then
              z1 = realone
              z2 = realone
              z3 = realone
              z4 = realone
           else
              z1 = realone
              z2 = -realone
              z3 = realone
              z4 = -realone
           end if
           lquery = lwork == -1
           if (m < 0) then
              info = -3
           else if (p < 0 .or. p > m) then
              info = -4
           else if (q < 0 .or. q > p .or. q > m - p .or. q > m - q) then
              info = -5
           else if (colmajor .and. ldx11 < max(1,p)) then
              info = -7
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
              info = -7
           else if (colmajor .and. ldx12 < max(1,p)) then
              info = -9
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
              info = -9
           else if (colmajor .and. ldx21 < max(1,m - p)) then
              info = -11
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
              info = -11
           else if (colmajor .and. ldx22 < max(1,m - p)) then
              info = -13
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
              info = -13
           end if
           ! compute workspace
           if (info == 0) then
              lworkopt = m - q
              lworkmin = m - q
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORBDB',-info)
              return
           else if (lquery) then
              return
           end if
           ! handle column-major and row-major separately
           if (colmajor) then
              ! reduce columns 1, ..., q of x11, x12, x21, and x22
              do i = 1,q
                 if (i == 1) then
                    call la_dscal(p - i + 1,z1,x11(i,i),1)
                 else
                    call la_dscal(p - i + 1,z1*cos(phi(i - 1)),x11(i,i),1)
                    call la_daxpy(p - i + 1,-z1*z3*z4*sin(phi(i - 1)),x12(i,i - 1),1,x11(i,i),1)

                 end if
                 if (i == 1) then
                    call la_dscal(m - p - i + 1,z2,x21(i,i),1)
                 else
                    call la_dscal(m - p - i + 1,z2*cos(phi(i - 1)),x21(i,i),1)
                    call la_daxpy(m - p - i + 1,-z2*z3*z4*sin(phi(i - 1)),x22(i,i - 1),1,x21(i,i), &
                              1)
                 end if
                 theta(i) = atan2(la_dnrm2(m - p - i + 1,x21(i,i),1),la_dnrm2(p - i + 1,x11( &
                           i,i),1))
                 if (p > i) then
                    call la_dlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
                 else if (p == i) then
                    call la_dlarfgp(p - i + 1,x11(i,i),x11(i,i),1,taup1(i))
                 end if
                 x11(i,i) = one
                 if (m - p > i) then
                    call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
                 else if (m - p == i) then
                    call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i,i),1,taup2(i))
                 end if
                 x21(i,i) = one
                 if (q > i) then
                    call la_dlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11, &
                              work)
                 end if
                 if (m - q + 1 > i) then
                    call la_dlarf('L',p - i + 1,m - q - i + 1,x11(i,i),1,taup1(i),x12(i,i),ldx12, &
                               work)
                 end if
                 if (q > i) then
                    call la_dlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21, &
                               work)
                 end if
                 if (m - q + 1 > i) then
                    call la_dlarf('L',m - p - i + 1,m - q - i + 1,x21(i,i),1,taup2(i),x22(i,i), &
                              ldx22,work)
                 end if
                 if (i < q) then
                    call la_dscal(q - i,-z1*z3*sin(theta(i)),x11(i,i + 1),ldx11)
                    call la_daxpy(q - i,z2*z3*cos(theta(i)),x21(i,i + 1),ldx21,x11(i,i + 1), &
                              ldx11)
                 end if
                 call la_dscal(m - q - i + 1,-z1*z4*sin(theta(i)),x12(i,i),ldx12)
                 call la_daxpy(m - q - i + 1,z2*z4*cos(theta(i)),x22(i,i),ldx22,x12(i,i),ldx12 &
                           )
                 if (i < q) phi(i) = atan2(la_dnrm2(q - i,x11(i,i + 1),ldx11),la_dnrm2( &
                           m - q - i + 1,x12(i,i),ldx12))
                 if (i < q) then
                    if (q - i == 1) then
                       call la_dlarfgp(q - i,x11(i,i + 1),x11(i,i + 1),ldx11,tauq1(i))
                    else
                       call la_dlarfgp(q - i,x11(i,i + 1),x11(i,i + 2),ldx11,tauq1(i))
                    end if
                    x11(i,i + 1) = one
                 end if
                 if (q + i - 1 < m) then
                    if (m - q == i) then
                       call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                    else
                       call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))

                    end if
                 end if
                 x12(i,i) = one
                 if (i < q) then
                    call la_dlarf('R',p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x11(i + 1,i + 1), &
                              ldx11,work)
                    call la_dlarf('R',m - p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 if (p > i) then
                    call la_dlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p > i) then
                    call la_dlarf('R',m - p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x22(i + 1,i), &
                              ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_dscal(m - q - i + 1,-z1*z4,x12(i,i),ldx12)
                 if (i >= m - q) then
                    call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                 else
                    call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))
                 end if
                 x12(i,i) = one
                 if (p > i) then
                    call la_dlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_dlarf('R',m - p - q,m - q - i + 1,x12(i,i),ldx12,tauq2(i), &
                            x22(q + 1,i),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_dscal(m - p - q - i + 1,z2*z4,x22(q + i,p + i),ldx22)
                 if (i == m - p - q) then
                    call la_dlarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i),ldx22,tauq2(p + i))

                 else
                    call la_dlarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i + 1),ldx22,tauq2(p + i) &
                               )
                 end if
                 x22(q + i,p + i) = one
                 if (i < m - p - q) then
                    call la_dlarf('R',m - p - q - i,m - p - q - i + 1,x22(q + i,p + i),ldx22,tauq2(p + i), &
                              x22(q + i + 1,p + i),ldx22,work)
                 end if
              end do
           else
              ! reduce columns 1, ..., q of x11, x12, x21, x22
              do i = 1,q
                 if (i == 1) then
                    call la_dscal(p - i + 1,z1,x11(i,i),ldx11)
                 else
                    call la_dscal(p - i + 1,z1*cos(phi(i - 1)),x11(i,i),ldx11)
                    call la_daxpy(p - i + 1,-z1*z3*z4*sin(phi(i - 1)),x12(i - 1,i),ldx12,x11(i,i), &
                               ldx11)
                 end if
                 if (i == 1) then
                    call la_dscal(m - p - i + 1,z2,x21(i,i),ldx21)
                 else
                    call la_dscal(m - p - i + 1,z2*cos(phi(i - 1)),x21(i,i),ldx21)
                    call la_daxpy(m - p - i + 1,-z2*z3*z4*sin(phi(i - 1)),x22(i - 1,i),ldx22,x21(i, &
                              i),ldx21)
                 end if
                 theta(i) = atan2(la_dnrm2(m - p - i + 1,x21(i,i),ldx21),la_dnrm2(p - i + 1, &
                           x11(i,i),ldx11))
                 call la_dlarfgp(p - i + 1,x11(i,i),x11(i,i + 1),ldx11,taup1(i))
                 x11(i,i) = one
                 if (i == m - p) then
                    call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i,i),ldx21,taup2(i))
                 else
                    call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i,i + 1),ldx21,taup2(i))
                 end if
                 x21(i,i) = one
                 if (q > i) then
                    call la_dlarf('R',q - i,p - i + 1,x11(i,i),ldx11,taup1(i),x11(i + 1,i), &
                              ldx11,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_dlarf('R',m - q - i + 1,p - i + 1,x11(i,i),ldx11,taup1(i),x12(i,i), &
                              ldx12,work)
                 end if
                 if (q > i) then
                    call la_dlarf('R',q - i,m - p - i + 1,x21(i,i),ldx21,taup2(i),x21(i + 1,i), &
                              ldx21,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_dlarf('R',m - q - i + 1,m - p - i + 1,x21(i,i),ldx21,taup2(i),x22(i,i), &
                              ldx22,work)
                 end if
                 if (i < q) then
                    call la_dscal(q - i,-z1*z3*sin(theta(i)),x11(i + 1,i),1)
                    call la_daxpy(q - i,z2*z3*cos(theta(i)),x21(i + 1,i),1,x11(i + 1,i),1)

                 end if
                 call la_dscal(m - q - i + 1,-z1*z4*sin(theta(i)),x12(i,i),1)
                 call la_daxpy(m - q - i + 1,z2*z4*cos(theta(i)),x22(i,i),1,x12(i,i),1)

                 if (i < q) phi(i) = atan2(la_dnrm2(q - i,x11(i + 1,i),1),la_dnrm2(m - q - &
                           i + 1,x12(i,i),1))
                 if (i < q) then
                    if (q - i == 1) then
                       call la_dlarfgp(q - i,x11(i + 1,i),x11(i + 1,i),1,tauq1(i))
                    else
                       call la_dlarfgp(q - i,x11(i + 1,i),x11(i + 2,i),1,tauq1(i))
                    end if
                    x11(i + 1,i) = one
                 end if
                 if (m - q > i) then
                    call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 else
                    call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i,i),1,tauq2(i))
                 end if
                 x12(i,i) = one
                 if (i < q) then
                    call la_dlarf('L',q - i,p - i,x11(i + 1,i),1,tauq1(i),x11(i + 1,i + 1),ldx11, &
                               work)
                    call la_dlarf('L',q - i,m - p - i,x11(i + 1,i),1,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 call la_dlarf('L',m - q - i + 1,p - i,x12(i,i),1,tauq2(i),x12(i,i + 1),ldx12, &
                           work)
                 if (m - p - i > 0) then
                    call la_dlarf('L',m - q - i + 1,m - p - i,x12(i,i),1,tauq2(i),x22(i,i + 1), &
                              ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_dscal(m - q - i + 1,-z1*z4,x12(i,i),1)
                 call la_dlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = one
                 if (p > i) then
                    call la_dlarf('L',m - q - i + 1,p - i,x12(i,i),1,tauq2(i),x12(i,i + 1),ldx12, &
                               work)
                 end if
                 if (m - p - q >= 1) call la_dlarf('L',m - q - i + 1,m - p - q,x12(i,i),1,tauq2(i), &
                           x22(i,q + 1),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_dscal(m - p - q - i + 1,z2*z4,x22(p + i,q + i),1)
                 if (m - p - q == i) then
                    call la_dlarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i,q + i),1,tauq2(p + i))

                 else
                    call la_dlarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i + 1,q + i),1,tauq2(p + i))

                    call la_dlarf('L',m - p - q - i + 1,m - p - q - i,x22(p + i,q + i),1,tauq2(p + i),x22(p + &
                              i,q + i + 1),ldx22,work)
                 end if
                 x22(p + i,q + i) = one
              end do
           end if
           return
     end subroutine la_dorbdb
     !> QORBDB: simultaneously bidiagonalizes the blocks of an M-by-M
     !> partitioned orthogonal matrix X:
     !> [ B11 | B12 0  0 ]
     !> [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T
     !> X = [-----------] = [---------] [----------------] [---------]   .
     !> [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
     !> [  0  |  0  0  I ]
     !> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     !> not the case, then X must be transposed and/or permuted. This can be
     !> done in constant time using the TRANS and SIGNS options. See QORCSD
     !> for details.)
     !> The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     !> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     !> represented implicitly by Householder vectors.
     !> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_qorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
               ldx22,theta,phi,taup1,taup2,tauq1,tauq2,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldx11,ldx12,ldx21,ldx22,lwork,m,p,q
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           real(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),tauq2(*),work(*)
           real(qp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ====================================================================
           ! Parameters
           real(qp),parameter :: realone = 1.0_qp

           ! Local Scalars
           logical(lk) :: colmajor,lquery
           integer(ilp) :: i,lworkmin,lworkopt
           real(qp) :: z1,z2,z3,z4
           ! Intrinsic Functions
           intrinsic :: atan2,cos,max,sin
           ! Executable Statements
           ! test input arguments
           info = 0
           colmajor = .not. la_lsame(trans,'T')
           if (.not. la_lsame(signs,'O')) then
              z1 = realone
              z2 = realone
              z3 = realone
              z4 = realone
           else
              z1 = realone
              z2 = -realone
              z3 = realone
              z4 = -realone
           end if
           lquery = lwork == -1
           if (m < 0) then
              info = -3
           else if (p < 0 .or. p > m) then
              info = -4
           else if (q < 0 .or. q > p .or. q > m - p .or. q > m - q) then
              info = -5
           else if (colmajor .and. ldx11 < max(1,p)) then
              info = -7
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
              info = -7
           else if (colmajor .and. ldx12 < max(1,p)) then
              info = -9
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
              info = -9
           else if (colmajor .and. ldx21 < max(1,m - p)) then
              info = -11
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
              info = -11
           else if (colmajor .and. ldx22 < max(1,m - p)) then
              info = -13
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
              info = -13
           end if
           ! compute workspace
           if (info == 0) then
              lworkopt = m - q
              lworkmin = m - q
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORBDB',-info)
              return
           else if (lquery) then
              return
           end if
           ! handle column-major and row-major separately
           if (colmajor) then
              ! reduce columns 1, ..., q of x11, x12, x21, and x22
              do i = 1,q
                 if (i == 1) then
                    call la_qscal(p - i + 1,z1,x11(i,i),1)
                 else
                    call la_qscal(p - i + 1,z1*cos(phi(i - 1)),x11(i,i),1)
                    call la_qaxpy(p - i + 1,-z1*z3*z4*sin(phi(i - 1)),x12(i,i - 1),1,x11(i,i),1)

                 end if
                 if (i == 1) then
                    call la_qscal(m - p - i + 1,z2,x21(i,i),1)
                 else
                    call la_qscal(m - p - i + 1,z2*cos(phi(i - 1)),x21(i,i),1)
                    call la_qaxpy(m - p - i + 1,-z2*z3*z4*sin(phi(i - 1)),x22(i,i - 1),1,x21(i,i), &
                              1)
                 end if
                 theta(i) = atan2(la_qnrm2(m - p - i + 1,x21(i,i),1),la_qnrm2(p - i + 1,x11( &
                           i,i),1))
                 if (p > i) then
                    call la_qlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
                 else if (p == i) then
                    call la_qlarfgp(p - i + 1,x11(i,i),x11(i,i),1,taup1(i))
                 end if
                 x11(i,i) = one
                 if (m - p > i) then
                    call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
                 else if (m - p == i) then
                    call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i,i),1,taup2(i))
                 end if
                 x21(i,i) = one
                 if (q > i) then
                    call la_qlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11, &
                              work)
                 end if
                 if (m - q + 1 > i) then
                    call la_qlarf('L',p - i + 1,m - q - i + 1,x11(i,i),1,taup1(i),x12(i,i),ldx12, &
                               work)
                 end if
                 if (q > i) then
                    call la_qlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21, &
                               work)
                 end if
                 if (m - q + 1 > i) then
                    call la_qlarf('L',m - p - i + 1,m - q - i + 1,x21(i,i),1,taup2(i),x22(i,i), &
                              ldx22,work)
                 end if
                 if (i < q) then
                    call la_qscal(q - i,-z1*z3*sin(theta(i)),x11(i,i + 1),ldx11)
                    call la_qaxpy(q - i,z2*z3*cos(theta(i)),x21(i,i + 1),ldx21,x11(i,i + 1), &
                              ldx11)
                 end if
                 call la_qscal(m - q - i + 1,-z1*z4*sin(theta(i)),x12(i,i),ldx12)
                 call la_qaxpy(m - q - i + 1,z2*z4*cos(theta(i)),x22(i,i),ldx22,x12(i,i),ldx12 &
                           )
                 if (i < q) phi(i) = atan2(la_qnrm2(q - i,x11(i,i + 1),ldx11),la_qnrm2( &
                           m - q - i + 1,x12(i,i),ldx12))
                 if (i < q) then
                    if (q - i == 1) then
                       call la_qlarfgp(q - i,x11(i,i + 1),x11(i,i + 1),ldx11,tauq1(i))
                    else
                       call la_qlarfgp(q - i,x11(i,i + 1),x11(i,i + 2),ldx11,tauq1(i))
                    end if
                    x11(i,i + 1) = one
                 end if
                 if (q + i - 1 < m) then
                    if (m - q == i) then
                       call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                    else
                       call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))

                    end if
                 end if
                 x12(i,i) = one
                 if (i < q) then
                    call la_qlarf('R',p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x11(i + 1,i + 1), &
                              ldx11,work)
                    call la_qlarf('R',m - p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 if (p > i) then
                    call la_qlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p > i) then
                    call la_qlarf('R',m - p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x22(i + 1,i), &
                              ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_qscal(m - q - i + 1,-z1*z4,x12(i,i),ldx12)
                 if (i >= m - q) then
                    call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                 else
                    call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))
                 end if
                 x12(i,i) = one
                 if (p > i) then
                    call la_qlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_qlarf('R',m - p - q,m - q - i + 1,x12(i,i),ldx12,tauq2(i), &
                            x22(q + 1,i),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_qscal(m - p - q - i + 1,z2*z4,x22(q + i,p + i),ldx22)
                 if (i == m - p - q) then
                    call la_qlarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i),ldx22,tauq2(p + i))

                 else
                    call la_qlarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i + 1),ldx22,tauq2(p + i) &
                               )
                 end if
                 x22(q + i,p + i) = one
                 if (i < m - p - q) then
                    call la_qlarf('R',m - p - q - i,m - p - q - i + 1,x22(q + i,p + i),ldx22,tauq2(p + i), &
                              x22(q + i + 1,p + i),ldx22,work)
                 end if
              end do
           else
              ! reduce columns 1, ..., q of x11, x12, x21, x22
              do i = 1,q
                 if (i == 1) then
                    call la_qscal(p - i + 1,z1,x11(i,i),ldx11)
                 else
                    call la_qscal(p - i + 1,z1*cos(phi(i - 1)),x11(i,i),ldx11)
                    call la_qaxpy(p - i + 1,-z1*z3*z4*sin(phi(i - 1)),x12(i - 1,i),ldx12,x11(i,i), &
                               ldx11)
                 end if
                 if (i == 1) then
                    call la_qscal(m - p - i + 1,z2,x21(i,i),ldx21)
                 else
                    call la_qscal(m - p - i + 1,z2*cos(phi(i - 1)),x21(i,i),ldx21)
                    call la_qaxpy(m - p - i + 1,-z2*z3*z4*sin(phi(i - 1)),x22(i - 1,i),ldx22,x21(i, &
                              i),ldx21)
                 end if
                 theta(i) = atan2(la_qnrm2(m - p - i + 1,x21(i,i),ldx21),la_qnrm2(p - i + 1, &
                           x11(i,i),ldx11))
                 call la_qlarfgp(p - i + 1,x11(i,i),x11(i,i + 1),ldx11,taup1(i))
                 x11(i,i) = one
                 if (i == m - p) then
                    call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i,i),ldx21,taup2(i))
                 else
                    call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i,i + 1),ldx21,taup2(i))
                 end if
                 x21(i,i) = one
                 if (q > i) then
                    call la_qlarf('R',q - i,p - i + 1,x11(i,i),ldx11,taup1(i),x11(i + 1,i), &
                              ldx11,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_qlarf('R',m - q - i + 1,p - i + 1,x11(i,i),ldx11,taup1(i),x12(i,i), &
                              ldx12,work)
                 end if
                 if (q > i) then
                    call la_qlarf('R',q - i,m - p - i + 1,x21(i,i),ldx21,taup2(i),x21(i + 1,i), &
                              ldx21,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_qlarf('R',m - q - i + 1,m - p - i + 1,x21(i,i),ldx21,taup2(i),x22(i,i), &
                              ldx22,work)
                 end if
                 if (i < q) then
                    call la_qscal(q - i,-z1*z3*sin(theta(i)),x11(i + 1,i),1)
                    call la_qaxpy(q - i,z2*z3*cos(theta(i)),x21(i + 1,i),1,x11(i + 1,i),1)

                 end if
                 call la_qscal(m - q - i + 1,-z1*z4*sin(theta(i)),x12(i,i),1)
                 call la_qaxpy(m - q - i + 1,z2*z4*cos(theta(i)),x22(i,i),1,x12(i,i),1)

                 if (i < q) phi(i) = atan2(la_qnrm2(q - i,x11(i + 1,i),1),la_qnrm2(m - q - &
                           i + 1,x12(i,i),1))
                 if (i < q) then
                    if (q - i == 1) then
                       call la_qlarfgp(q - i,x11(i + 1,i),x11(i + 1,i),1,tauq1(i))
                    else
                       call la_qlarfgp(q - i,x11(i + 1,i),x11(i + 2,i),1,tauq1(i))
                    end if
                    x11(i + 1,i) = one
                 end if
                 if (m - q > i) then
                    call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 else
                    call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i,i),1,tauq2(i))
                 end if
                 x12(i,i) = one
                 if (i < q) then
                    call la_qlarf('L',q - i,p - i,x11(i + 1,i),1,tauq1(i),x11(i + 1,i + 1),ldx11, &
                               work)
                    call la_qlarf('L',q - i,m - p - i,x11(i + 1,i),1,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 call la_qlarf('L',m - q - i + 1,p - i,x12(i,i),1,tauq2(i),x12(i,i + 1),ldx12, &
                           work)
                 if (m - p - i > 0) then
                    call la_qlarf('L',m - q - i + 1,m - p - i,x12(i,i),1,tauq2(i),x22(i,i + 1), &
                              ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_qscal(m - q - i + 1,-z1*z4,x12(i,i),1)
                 call la_qlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = one
                 if (p > i) then
                    call la_qlarf('L',m - q - i + 1,p - i,x12(i,i),1,tauq2(i),x12(i,i + 1),ldx12, &
                               work)
                 end if
                 if (m - p - q >= 1) call la_qlarf('L',m - q - i + 1,m - p - q,x12(i,i),1,tauq2(i), &
                           x22(i,q + 1),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_qscal(m - p - q - i + 1,z2*z4,x22(p + i,q + i),1)
                 if (m - p - q == i) then
                    call la_qlarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i,q + i),1,tauq2(p + i))

                 else
                    call la_qlarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i + 1,q + i),1,tauq2(p + i))

                    call la_qlarf('L',m - p - q - i + 1,m - p - q - i,x22(p + i,q + i),1,tauq2(p + i),x22(p + &
                              i,q + i + 1),ldx22,work)
                 end if
                 x22(p + i,q + i) = one
              end do
           end if
           return
     end subroutine la_qorbdb

     !> SORBDB5: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then some other vector from the orthogonal complement
     !> is returned. This vector is chosen in an arbitrary but deterministic
     !> way.

     pure subroutine la_sorbdb5(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           real(sp),intent(out) :: work(*)
           real(sp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,j
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('SORBDB5',-info)
              return
           end if
           ! project x onto the orthogonal complement of q
           call la_sorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work,lwork, &
                     childinfo)
           ! if the projection is nonzero, then return
           if (la_snrm2(m1,x1,incx1) /= zero .or. la_snrm2(m2,x2,incx2) /= zero) &
                     then
              return
           end if
           ! project each standard basis vector e_1,...,e_m1 in turn, stopping
           ! when a nonzero projection is found
           do i = 1,m1
              do j = 1,m1
                 x1(j) = zero
              end do
              x1(i) = one
              do j = 1,m2
                 x2(j) = zero
              end do
              call la_sorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_snrm2(m1,x1,incx1) /= zero .or. la_snrm2(m2,x2,incx2) /= zero) &
                        then
                 return
              end if
           end do
           ! project each standard basis vector e_(m1+1),...,e_(m1+m2) in turn,
           ! stopping when a nonzero projection is found
           do i = 1,m2
              do j = 1,m1
                 x1(j) = zero
              end do
              do j = 1,m2
                 x2(j) = zero
              end do
              x2(i) = one
              call la_sorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_snrm2(m1,x1,incx1) /= zero .or. la_snrm2(m2,x2,incx2) /= zero) &
                        then
                 return
              end if
           end do
           return
     end subroutine la_sorbdb5
     !> DORBDB5: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then some other vector from the orthogonal complement
     !> is returned. This vector is chosen in an arbitrary but deterministic
     !> way.

     pure subroutine la_dorbdb5(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           real(dp),intent(out) :: work(*)
           real(dp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,j
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('DORBDB5',-info)
              return
           end if
           ! project x onto the orthogonal complement of q
           call la_dorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work,lwork, &
                     childinfo)
           ! if the projection is nonzero, then return
           if (la_dnrm2(m1,x1,incx1) /= zero .or. la_dnrm2(m2,x2,incx2) /= zero) &
                     then
              return
           end if
           ! project each standard basis vector e_1,...,e_m1 in turn, stopping
           ! when a nonzero projection is found
           do i = 1,m1
              do j = 1,m1
                 x1(j) = zero
              end do
              x1(i) = one
              do j = 1,m2
                 x2(j) = zero
              end do
              call la_dorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_dnrm2(m1,x1,incx1) /= zero .or. la_dnrm2(m2,x2,incx2) /= zero) &
                        then
                 return
              end if
           end do
           ! project each standard basis vector e_(m1+1),...,e_(m1+m2) in turn,
           ! stopping when a nonzero projection is found
           do i = 1,m2
              do j = 1,m1
                 x1(j) = zero
              end do
              do j = 1,m2
                 x2(j) = zero
              end do
              x2(i) = one
              call la_dorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_dnrm2(m1,x1,incx1) /= zero .or. la_dnrm2(m2,x2,incx2) /= zero) &
                        then
                 return
              end if
           end do
           return
     end subroutine la_dorbdb5
     !> QORBDB5: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then some other vector from the orthogonal complement
     !> is returned. This vector is chosen in an arbitrary but deterministic
     !> way.

     pure subroutine la_qorbdb5(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           real(qp),intent(out) :: work(*)
           real(qp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,j
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('QORBDB5',-info)
              return
           end if
           ! project x onto the orthogonal complement of q
           call la_qorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work,lwork, &
                     childinfo)
           ! if the projection is nonzero, then return
           if (la_qnrm2(m1,x1,incx1) /= zero .or. la_qnrm2(m2,x2,incx2) /= zero) &
                     then
              return
           end if
           ! project each standard basis vector e_1,...,e_m1 in turn, stopping
           ! when a nonzero projection is found
           do i = 1,m1
              do j = 1,m1
                 x1(j) = zero
              end do
              x1(i) = one
              do j = 1,m2
                 x2(j) = zero
              end do
              call la_qorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_qnrm2(m1,x1,incx1) /= zero .or. la_qnrm2(m2,x2,incx2) /= zero) &
                        then
                 return
              end if
           end do
           ! project each standard basis vector e_(m1+1),...,e_(m1+m2) in turn,
           ! stopping when a nonzero projection is found
           do i = 1,m2
              do j = 1,m1
                 x1(j) = zero
              end do
              do j = 1,m2
                 x2(j) = zero
              end do
              x2(i) = one
              call la_qorbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_qnrm2(m1,x1,incx1) /= zero .or. la_qnrm2(m2,x2,incx2) /= zero) &
                        then
                 return
              end if
           end do
           return
     end subroutine la_qorbdb5

     !> SORCSD: computes the CS decomposition of an M-by-M partitioned
     !> orthogonal matrix X:
     !> [  I  0  0 |  0  0  0 ]
     !> [  0  C  0 |  0 -S  0 ]
     !> [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T
     !> X = [-----------] = [---------] [---------------------] [---------]   .
     !> [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
     !> [  0  S  0 |  0  C  0 ]
     !> [  0  0  I |  0  0  0 ]
     !> X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P,
     !> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     !> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     !> which R = MIN(P,M-P,Q,M-Q).

     recursive subroutine la_sorcsd(jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,p,q,x11, &
     ldx11,x12,ldx12,x21,ldx21,x22,ldx22,theta,u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t, &
               work,lwork,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,ldx11,ldx12,ldx21,ldx22, &
                     lwork,m,p,q
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: theta(*)
           real(sp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*),work(*)

           real(sp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ===================================================================

           ! Local Arrays
           real(sp) :: dummy(1)
           ! Local Scalars
           character :: transt,signst
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,itauq2,j,lbbcsdwork, &
           lbbcsdworkmin,lbbcsdworkopt,lorbdbwork,lorbdbworkmin,lorbdbworkopt,lorglqwork, &
           lorglqworkmin,lorglqworkopt,lorgqrwork,lorgqrworkmin,lorgqrworkopt,lworkmin, &
                     lworkopt
           logical(lk) :: colmajor,defaultsigns,lquery,wantu1,wantu2,wantv1t,wantv2t
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           defaultsigns = .not. la_lsame(signs,'O')
           lquery = lwork == -1
           if (m < 0) then
              info = -7
           else if (p < 0 .or. p > m) then
              info = -8
           else if (q < 0 .or. q > m) then
              info = -9
           else if (colmajor .and. ldx11 < max(1,p)) then
             info = -11
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
             info = -11
           else if (colmajor .and. ldx12 < max(1,p)) then
             info = -13
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
             info = -13
           else if (colmajor .and. ldx21 < max(1,m - p)) then
             info = -15
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
             info = -15
           else if (colmajor .and. ldx22 < max(1,m - p)) then
             info = -17
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
             info = -17
           else if (wantu1 .and. ldu1 < p) then
              info = -20
           else if (wantu2 .and. ldu2 < m - p) then
              info = -22
           else if (wantv1t .and. ldv1t < q) then
              info = -24
           else if (wantv2t .and. ldv2t < m - q) then
              info = -26
           end if
           ! work with transpose if convenient
           if (info == 0 .and. min(p,m - p) < min(q,m - q)) then
              if (colmajor) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_sorcsd(jobv1t,jobv2t,jobu1,jobu2,transt,signst,m,q,p,x11, &
              ldx11,x21,ldx21,x12,ldx12,x22,ldx22,theta,v1t,ldv1t,v2t,ldv2t,u1,ldu1, &
                        u2,ldu2,work,lwork,iwork,info)
              return
           end if
           ! work with permutation [ 0 i; i 0 ] * x * [ 0 i; i 0 ] if
           ! convenient
           if (info == 0 .and. m - q < q) then
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_sorcsd(jobu2,jobu1,jobv2t,jobv1t,trans,signst,m,m - p,m - q,x22, &
              ldx22,x21,ldx21,x12,ldx12,x11,ldx11,theta,u2,ldu2,u1,ldu1,v2t,ldv2t, &
                        v1t,ldv1t,work,lwork,iwork,info)
              return
           end if
           ! compute workspace
           if (info == 0) then
              iphi = 2
              itaup1 = iphi + max(1,q - 1)
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              itauq2 = itauq1 + max(1,q)
              iorgqr = itauq2 + max(1,m - q)
              call la_sorgqr(m - q,m - q,m - q,dummy,max(1,m - q),dummy,work,-1,childinfo)

              lorgqrworkopt = int(work(1),KIND=ilp)
              lorgqrworkmin = max(1,m - q)
              iorglq = itauq2 + max(1,m - q)
              call la_sorglq(m - q,m - q,m - q,dummy,max(1,m - q),dummy,work,-1,childinfo)

              lorglqworkopt = int(work(1),KIND=ilp)
              lorglqworkmin = max(1,m - q)
              iorbdb = itauq2 + max(1,m - q)
              call la_sorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
                        ldx22,dummy,dummy,dummy,dummy,dummy,dummy,work,-1,childinfo)
              lorbdbworkopt = int(work(1),KIND=ilp)
              lorbdbworkmin = lorbdbworkopt
              ib11d = itauq2 + max(1,m - q)
              ib11e = ib11d + max(1,q)
              ib12d = ib11e + max(1,q - 1)
              ib12e = ib12d + max(1,q)
              ib21d = ib12e + max(1,q - 1)
              ib21e = ib21d + max(1,q)
              ib22d = ib21e + max(1,q - 1)
              ib22e = ib22d + max(1,q)
              ibbcsd = ib22e + max(1,q - 1)
              call la_sbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,dummy,dummy,u1, &
              ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,dummy,dummy,dummy,dummy,dummy,dummy, &
                        dummy,dummy,work,-1,childinfo)
              lbbcsdworkopt = int(work(1),KIND=ilp)
              lbbcsdworkmin = lbbcsdworkopt
              lworkopt = max(iorgqr + lorgqrworkopt,iorglq + lorglqworkopt,iorbdb + &
                        lorbdbworkopt,ibbcsd + lbbcsdworkopt) - 1
              lworkmin = max(iorgqr + lorgqrworkmin,iorglq + lorglqworkmin,iorbdb + &
                        lorbdbworkopt,ibbcsd + lbbcsdworkmin) - 1
              work(1) = max(lworkopt,lworkmin)
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -22
              else
                 lorgqrwork = lwork - iorgqr + 1
                 lorglqwork = lwork - iorglq + 1
                 lorbdbwork = lwork - iorbdb + 1
                 lbbcsdwork = lwork - ibbcsd + 1
              end if
           end if
           ! abort if any illegal arguments
           if (info /= 0) then
              call la_xerbla('SORCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! transform to bidiagonal block form
           call la_sorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
           ldx22,theta,work(iphi),work(itaup1),work(itaup2),work(itauq1),work(itauq2),work( &
                     iorbdb),lorbdbwork,childinfo)
           ! accumulate householder reflectors
           if (colmajor) then
              if (wantu1 .and. p > 0) then
                 call la_slacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_sorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqrwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_slacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_sorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqrwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_slacpy('U',q - 1,q - 1,x11(1,2),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_sorglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglqwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_slacpy('U',p,m - q,x12,ldx12,v2t,ldv2t)
                 call la_slacpy('U',m - p - q,m - p - q,x22(q + 1,p + 1),ldx22,v2t(p + 1,p + 1),ldv2t)

                    call la_sorglq(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorglq), &
                              lorglqwork,info)
              end if
           else
              if (wantu1 .and. p > 0) then
                 call la_slacpy('U',q,p,x11,ldx11,u1,ldu1)
                 call la_sorglq(p,p,q,u1,ldu1,work(itaup1),work(iorglq),lorglqwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_slacpy('U',q,m - p,x21,ldx21,u2,ldu2)
                 call la_sorglq(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorglq),lorglqwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_slacpy('L',q - 1,q - 1,x11(2,1),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_sorgqr(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorgqr), &
                           lorgqrwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_slacpy('L',m - q,p,x12,ldx12,v2t,ldv2t)
                 call la_slacpy('L',m - p - q,m - p - q,x22(p + 1,q + 1),ldx22,v2t(p + 1,p + 1),ldv2t)

                 call la_sorgqr(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorgqr), &
                           lorgqrwork,info)
              end if
           end if
           ! compute the csd of the matrix in bidiagonal-block form
           call la_sbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,work(iphi),u1, &
            ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,work(ib11d),work(ib11e),work(ib12d),work( &
            ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsdwork, &
                      info)
           ! permute rows and columns to place identity submatrices in top-
           ! left corner of (1,1)-block and/or bottom-right corner of (1,2)-
           ! block and/or bottom-right corner of (2,1)-block and/or top-left
           ! corner of (2,2)-block
           if (q > 0 .and. wantu2) then
              do i = 1,q
                 iwork(i) = m - p - q + i
              end do
              do i = q + 1,m - p
                 iwork(i) = i - q
              end do
              if (colmajor) then
                 call la_slapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              else
                 call la_slapmr(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           end if
           if (m > 0 .and. wantv2t) then
              do i = 1,p
                 iwork(i) = m - p - q + i
              end do
              do i = p + 1,m - q
                 iwork(i) = i - p
              end do
              if (.not. colmajor) then
                 call la_slapmt(.false.,m - q,m - q,v2t,ldv2t,iwork)
              else
                 call la_slapmr(.false.,m - q,m - q,v2t,ldv2t,iwork)
              end if
           end if
           return
           ! end la_sorcsd
     end subroutine la_sorcsd
     !> DORCSD: computes the CS decomposition of an M-by-M partitioned
     !> orthogonal matrix X:
     !> [  I  0  0 |  0  0  0 ]
     !> [  0  C  0 |  0 -S  0 ]
     !> [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T
     !> X = [-----------] = [---------] [---------------------] [---------]   .
     !> [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
     !> [  0  S  0 |  0  C  0 ]
     !> [  0  0  I |  0  0  0 ]
     !> X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P,
     !> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     !> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     !> which R = MIN(P,M-P,Q,M-Q).

     recursive subroutine la_dorcsd(jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,p,q,x11, &
     ldx11,x12,ldx12,x21,ldx21,x22,ldx22,theta,u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t, &
               work,lwork,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,ldx11,ldx12,ldx21,ldx22, &
                     lwork,m,p,q
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: theta(*)
           real(dp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*),work(*)

           real(dp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ===================================================================

           ! Local Scalars
           character :: transt,signst
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,itauq2,j,lbbcsdwork, &
           lbbcsdworkmin,lbbcsdworkopt,lorbdbwork,lorbdbworkmin,lorbdbworkopt,lorglqwork, &
           lorglqworkmin,lorglqworkopt,lorgqrwork,lorgqrworkmin,lorgqrworkopt,lworkmin, &
                     lworkopt
           logical(lk) :: colmajor,defaultsigns,lquery,wantu1,wantu2,wantv1t,wantv2t
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           defaultsigns = .not. la_lsame(signs,'O')
           lquery = lwork == -1
           if (m < 0) then
              info = -7
           else if (p < 0 .or. p > m) then
              info = -8
           else if (q < 0 .or. q > m) then
              info = -9
           else if (colmajor .and. ldx11 < max(1,p)) then
             info = -11
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
             info = -11
           else if (colmajor .and. ldx12 < max(1,p)) then
             info = -13
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
             info = -13
           else if (colmajor .and. ldx21 < max(1,m - p)) then
             info = -15
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
             info = -15
           else if (colmajor .and. ldx22 < max(1,m - p)) then
             info = -17
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
             info = -17
           else if (wantu1 .and. ldu1 < p) then
              info = -20
           else if (wantu2 .and. ldu2 < m - p) then
              info = -22
           else if (wantv1t .and. ldv1t < q) then
              info = -24
           else if (wantv2t .and. ldv2t < m - q) then
              info = -26
           end if
           ! work with transpose if convenient
           if (info == 0 .and. min(p,m - p) < min(q,m - q)) then
              if (colmajor) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_dorcsd(jobv1t,jobv2t,jobu1,jobu2,transt,signst,m,q,p,x11, &
              ldx11,x21,ldx21,x12,ldx12,x22,ldx22,theta,v1t,ldv1t,v2t,ldv2t,u1,ldu1, &
                        u2,ldu2,work,lwork,iwork,info)
              return
           end if
           ! work with permutation [ 0 i; i 0 ] * x * [ 0 i; i 0 ] if
           ! convenient
           if (info == 0 .and. m - q < q) then
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_dorcsd(jobu2,jobu1,jobv2t,jobv1t,trans,signst,m,m - p,m - q,x22, &
              ldx22,x21,ldx21,x12,ldx12,x11,ldx11,theta,u2,ldu2,u1,ldu1,v2t,ldv2t, &
                        v1t,ldv1t,work,lwork,iwork,info)
              return
           end if
           ! compute workspace
           if (info == 0) then
              iphi = 2
              itaup1 = iphi + max(1,q - 1)
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              itauq2 = itauq1 + max(1,q)
              iorgqr = itauq2 + max(1,m - q)
              call la_dorgqr(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorgqrworkopt = int(work(1),KIND=ilp)
              lorgqrworkmin = max(1,m - q)
              iorglq = itauq2 + max(1,m - q)
              call la_dorglq(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorglqworkopt = int(work(1),KIND=ilp)
              lorglqworkmin = max(1,m - q)
              iorbdb = itauq2 + max(1,m - q)
              call la_dorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
                        ldx22,theta,v1t,u1,u2,v1t,v2t,work,-1,childinfo)
              lorbdbworkopt = int(work(1),KIND=ilp)
              lorbdbworkmin = lorbdbworkopt
              ib11d = itauq2 + max(1,m - q)
              ib11e = ib11d + max(1,q)
              ib12d = ib11e + max(1,q - 1)
              ib12e = ib12d + max(1,q)
              ib21d = ib12e + max(1,q - 1)
              ib21e = ib21d + max(1,q)
              ib22d = ib21e + max(1,q - 1)
              ib22e = ib22d + max(1,q)
              ibbcsd = ib22e + max(1,q - 1)
              call la_dbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,theta,u1, &
              ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,u1,u1,u1,u1,u1,u1,u1,u1,work,-1, &
                        childinfo)
              lbbcsdworkopt = int(work(1),KIND=ilp)
              lbbcsdworkmin = lbbcsdworkopt
              lworkopt = max(iorgqr + lorgqrworkopt,iorglq + lorglqworkopt,iorbdb + &
                        lorbdbworkopt,ibbcsd + lbbcsdworkopt) - 1
              lworkmin = max(iorgqr + lorgqrworkmin,iorglq + lorglqworkmin,iorbdb + &
                        lorbdbworkopt,ibbcsd + lbbcsdworkmin) - 1
              work(1) = max(lworkopt,lworkmin)
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -22
              else
                 lorgqrwork = lwork - iorgqr + 1
                 lorglqwork = lwork - iorglq + 1
                 lorbdbwork = lwork - iorbdb + 1
                 lbbcsdwork = lwork - ibbcsd + 1
              end if
           end if
           ! abort if any illegal arguments
           if (info /= 0) then
              call la_xerbla('DORCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! transform to bidiagonal block form
           call la_dorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
           ldx22,theta,work(iphi),work(itaup1),work(itaup2),work(itauq1),work(itauq2),work( &
                     iorbdb),lorbdbwork,childinfo)
           ! accumulate householder reflectors
           if (colmajor) then
              if (wantu1 .and. p > 0) then
                 call la_dlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_dorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqrwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_dlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_dorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqrwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_dlacpy('U',q - 1,q - 1,x11(1,2),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_dorglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglqwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_dlacpy('U',p,m - q,x12,ldx12,v2t,ldv2t)
                 if (m - p > q) then
                    call la_dlacpy('U',m - p - q,m - p - q,x22(q + 1,p + 1),ldx22,v2t(p + 1,p + 1), &
                              ldv2t)
                 end if
                 if (m > q) then
                    call la_dorglq(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorglq), &
                              lorglqwork,info)
                 end if
              end if
           else
              if (wantu1 .and. p > 0) then
                 call la_dlacpy('U',q,p,x11,ldx11,u1,ldu1)
                 call la_dorglq(p,p,q,u1,ldu1,work(itaup1),work(iorglq),lorglqwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_dlacpy('U',q,m - p,x21,ldx21,u2,ldu2)
                 call la_dorglq(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorglq),lorglqwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_dlacpy('L',q - 1,q - 1,x11(2,1),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_dorgqr(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorgqr), &
                           lorgqrwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_dlacpy('L',m - q,p,x12,ldx12,v2t,ldv2t)
                 call la_dlacpy('L',m - p - q,m - p - q,x22(p + 1,q + 1),ldx22,v2t(p + 1,p + 1),ldv2t)

                 call la_dorgqr(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorgqr), &
                           lorgqrwork,info)
              end if
           end if
           ! compute the csd of the matrix in bidiagonal-block form
           call la_dbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,work(iphi),u1, &
            ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,work(ib11d),work(ib11e),work(ib12d),work( &
            ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsdwork, &
                      info)
           ! permute rows and columns to place identity submatrices in top-
           ! left corner of (1,1)-block and/or bottom-right corner of (1,2)-
           ! block and/or bottom-right corner of (2,1)-block and/or top-left
           ! corner of (2,2)-block
           if (q > 0 .and. wantu2) then
              do i = 1,q
                 iwork(i) = m - p - q + i
              end do
              do i = q + 1,m - p
                 iwork(i) = i - q
              end do
              if (colmajor) then
                 call la_dlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              else
                 call la_dlapmr(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           end if
           if (m > 0 .and. wantv2t) then
              do i = 1,p
                 iwork(i) = m - p - q + i
              end do
              do i = p + 1,m - q
                 iwork(i) = i - p
              end do
              if (.not. colmajor) then
                 call la_dlapmt(.false.,m - q,m - q,v2t,ldv2t,iwork)
              else
                 call la_dlapmr(.false.,m - q,m - q,v2t,ldv2t,iwork)
              end if
           end if
           return
           ! end la_dorcsd
     end subroutine la_dorcsd
     !> QORCSD: computes the CS decomposition of an M-by-M partitioned
     !> orthogonal matrix X:
     !> [  I  0  0 |  0  0  0 ]
     !> [  0  C  0 |  0 -S  0 ]
     !> [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T
     !> X = [-----------] = [---------] [---------------------] [---------]   .
     !> [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
     !> [  0  S  0 |  0  C  0 ]
     !> [  0  0  I |  0  0  0 ]
     !> X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P,
     !> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     !> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     !> which R = MIN(P,M-P,Q,M-Q).

     recursive subroutine la_qorcsd(jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,p,q,x11, &
     ldx11,x12,ldx12,x21,ldx21,x22,ldx22,theta,u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t, &
               work,lwork,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,ldx11,ldx12,ldx21,ldx22, &
                     lwork,m,p,q
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: theta(*)
           real(qp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*),work(*)

           real(qp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ===================================================================

           ! Local Scalars
           character :: transt,signst
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,itauq2,j,lbbcsdwork, &
           lbbcsdworkmin,lbbcsdworkopt,lorbdbwork,lorbdbworkmin,lorbdbworkopt,lorglqwork, &
           lorglqworkmin,lorglqworkopt,lorgqrwork,lorgqrworkmin,lorgqrworkopt,lworkmin, &
                     lworkopt
           logical(lk) :: colmajor,defaultsigns,lquery,wantu1,wantu2,wantv1t,wantv2t
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           defaultsigns = .not. la_lsame(signs,'O')
           lquery = lwork == -1
           if (m < 0) then
              info = -7
           else if (p < 0 .or. p > m) then
              info = -8
           else if (q < 0 .or. q > m) then
              info = -9
           else if (colmajor .and. ldx11 < max(1,p)) then
             info = -11
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
             info = -11
           else if (colmajor .and. ldx12 < max(1,p)) then
             info = -13
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
             info = -13
           else if (colmajor .and. ldx21 < max(1,m - p)) then
             info = -15
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
             info = -15
           else if (colmajor .and. ldx22 < max(1,m - p)) then
             info = -17
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
             info = -17
           else if (wantu1 .and. ldu1 < p) then
              info = -20
           else if (wantu2 .and. ldu2 < m - p) then
              info = -22
           else if (wantv1t .and. ldv1t < q) then
              info = -24
           else if (wantv2t .and. ldv2t < m - q) then
              info = -26
           end if
           ! work with transpose if convenient
           if (info == 0 .and. min(p,m - p) < min(q,m - q)) then
              if (colmajor) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_qorcsd(jobv1t,jobv2t,jobu1,jobu2,transt,signst,m,q,p,x11, &
              ldx11,x21,ldx21,x12,ldx12,x22,ldx22,theta,v1t,ldv1t,v2t,ldv2t,u1,ldu1, &
                        u2,ldu2,work,lwork,iwork,info)
              return
           end if
           ! work with permutation [ 0 i; i 0 ] * x * [ 0 i; i 0 ] if
           ! convenient
           if (info == 0 .and. m - q < q) then
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_qorcsd(jobu2,jobu1,jobv2t,jobv1t,trans,signst,m,m - p,m - q,x22, &
              ldx22,x21,ldx21,x12,ldx12,x11,ldx11,theta,u2,ldu2,u1,ldu1,v2t,ldv2t, &
                        v1t,ldv1t,work,lwork,iwork,info)
              return
           end if
           ! compute workspace
           if (info == 0) then
              iphi = 2
              itaup1 = iphi + max(1,q - 1)
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              itauq2 = itauq1 + max(1,q)
              iorgqr = itauq2 + max(1,m - q)
              call la_qorgqr(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorgqrworkopt = int(work(1),KIND=ilp)
              lorgqrworkmin = max(1,m - q)
              iorglq = itauq2 + max(1,m - q)
              call la_qorglq(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorglqworkopt = int(work(1),KIND=ilp)
              lorglqworkmin = max(1,m - q)
              iorbdb = itauq2 + max(1,m - q)
              call la_qorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
                        ldx22,theta,v1t,u1,u2,v1t,v2t,work,-1,childinfo)
              lorbdbworkopt = int(work(1),KIND=ilp)
              lorbdbworkmin = lorbdbworkopt
              ib11d = itauq2 + max(1,m - q)
              ib11e = ib11d + max(1,q)
              ib12d = ib11e + max(1,q - 1)
              ib12e = ib12d + max(1,q)
              ib21d = ib12e + max(1,q - 1)
              ib21e = ib21d + max(1,q)
              ib22d = ib21e + max(1,q - 1)
              ib22e = ib22d + max(1,q)
              ibbcsd = ib22e + max(1,q - 1)
              call la_qbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,theta,u1, &
              ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,u1,u1,u1,u1,u1,u1,u1,u1,work,-1, &
                        childinfo)
              lbbcsdworkopt = int(work(1),KIND=ilp)
              lbbcsdworkmin = lbbcsdworkopt
              lworkopt = max(iorgqr + lorgqrworkopt,iorglq + lorglqworkopt,iorbdb + &
                        lorbdbworkopt,ibbcsd + lbbcsdworkopt) - 1
              lworkmin = max(iorgqr + lorgqrworkmin,iorglq + lorglqworkmin,iorbdb + &
                        lorbdbworkopt,ibbcsd + lbbcsdworkmin) - 1
              work(1) = max(lworkopt,lworkmin)
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -22
              else
                 lorgqrwork = lwork - iorgqr + 1
                 lorglqwork = lwork - iorglq + 1
                 lorbdbwork = lwork - iorbdb + 1
                 lbbcsdwork = lwork - ibbcsd + 1
              end if
           end if
           ! abort if any illegal arguments
           if (info /= 0) then
              call la_xerbla('QORCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! transform to bidiagonal block form
           call la_qorbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
           ldx22,theta,work(iphi),work(itaup1),work(itaup2),work(itauq1),work(itauq2),work( &
                     iorbdb),lorbdbwork,childinfo)
           ! accumulate householder reflectors
           if (colmajor) then
              if (wantu1 .and. p > 0) then
                 call la_qlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_qorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqrwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_qlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_qorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqrwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_qlacpy('U',q - 1,q - 1,x11(1,2),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_qorglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglqwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_qlacpy('U',p,m - q,x12,ldx12,v2t,ldv2t)
                 if (m - p > q) then
                    call la_qlacpy('U',m - p - q,m - p - q,x22(q + 1,p + 1),ldx22,v2t(p + 1,p + 1), &
                              ldv2t)
                 end if
                 if (m > q) then
                    call la_qorglq(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorglq), &
                              lorglqwork,info)
                 end if
              end if
           else
              if (wantu1 .and. p > 0) then
                 call la_qlacpy('U',q,p,x11,ldx11,u1,ldu1)
                 call la_qorglq(p,p,q,u1,ldu1,work(itaup1),work(iorglq),lorglqwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_qlacpy('U',q,m - p,x21,ldx21,u2,ldu2)
                 call la_qorglq(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorglq),lorglqwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_qlacpy('L',q - 1,q - 1,x11(2,1),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_qorgqr(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorgqr), &
                           lorgqrwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_qlacpy('L',m - q,p,x12,ldx12,v2t,ldv2t)
                 call la_qlacpy('L',m - p - q,m - p - q,x22(p + 1,q + 1),ldx22,v2t(p + 1,p + 1),ldv2t)

                 call la_qorgqr(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorgqr), &
                           lorgqrwork,info)
              end if
           end if
           ! compute the csd of the matrix in bidiagonal-block form
           call la_qbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,work(iphi),u1, &
            ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,work(ib11d),work(ib11e),work(ib12d),work( &
            ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsdwork, &
                      info)
           ! permute rows and columns to place identity submatrices in top-
           ! left corner of (1,1)-block and/or bottom-right corner of (1,2)-
           ! block and/or bottom-right corner of (2,1)-block and/or top-left
           ! corner of (2,2)-block
           if (q > 0 .and. wantu2) then
              do i = 1,q
                 iwork(i) = m - p - q + i
              end do
              do i = q + 1,m - p
                 iwork(i) = i - q
              end do
              if (colmajor) then
                 call la_qlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              else
                 call la_qlapmr(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           end if
           if (m > 0 .and. wantv2t) then
              do i = 1,p
                 iwork(i) = m - p - q + i
              end do
              do i = p + 1,m - q
                 iwork(i) = i - p
              end do
              if (.not. colmajor) then
                 call la_qlapmt(.false.,m - q,m - q,v2t,ldv2t,iwork)
              else
                 call la_qlapmr(.false.,m - q,m - q,v2t,ldv2t,iwork)
              end if
           end if
           return
           ! end la_qorcsd
     end subroutine la_qorcsd

     !> SORBDB1: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
     !> M-P, or M-Q. Routines SORBDB2, SORBDB3, and SORBDB4 handle cases in
     !> which Q is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_sorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           real(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < q .or. m - p < q) then
              info = -2
           else if (q < 0 .or. m - q < q) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 2
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SORBDB1',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., q of x11 and x21
           do i = 1,q
              call la_slarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              call la_slarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              theta(i) = atan2(x21(i,i),x11(i,i))
              c = cos(theta(i))
              s = sin(theta(i))
              x11(i,i) = one
              x21(i,i) = one
              call la_slarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
              call la_slarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
              if (i < q) then
                 call la_srot(q - i,x11(i,i + 1),ldx11,x21(i,i + 1),ldx21,c,s)
                 call la_slarfgp(q - i,x21(i,i + 1),x21(i,i + 2),ldx21,tauq1(i))
                 s = x21(i,i + 1)
                 x21(i,i + 1) = one
                 call la_slarf('R',p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
                 call la_slarf('R',m - p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
                 c = sqrt(la_snrm2(p - i,x11(i + 1,i + 1),1)**2 + la_snrm2(m - p - i,x21(i + 1, &
                           i + 1),1)**2)
                 phi(i) = atan2(s,c)
                 call la_sorbdb5(p - i,m - p - i,q - i - 1,x11(i + 1,i + 1),1,x21(i + 1,i + 1),1,x11(i + 1, &
                           i + 2),ldx11,x21(i + 1,i + 2),ldx21,work(iorbdb5),lorbdb5,childinfo)
              end if
           end do
           return
     end subroutine la_sorbdb1
     !> DORBDB1: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
     !> M-P, or M-Q. Routines DORBDB2, DORBDB3, and DORBDB4 handle cases in
     !> which Q is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_dorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           real(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < q .or. m - p < q) then
              info = -2
           else if (q < 0 .or. m - q < q) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 2
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DORBDB1',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., q of x11 and x21
           do i = 1,q
              call la_dlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              theta(i) = atan2(x21(i,i),x11(i,i))
              c = cos(theta(i))
              s = sin(theta(i))
              x11(i,i) = one
              x21(i,i) = one
              call la_dlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
              call la_dlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
              if (i < q) then
                 call la_drot(q - i,x11(i,i + 1),ldx11,x21(i,i + 1),ldx21,c,s)
                 call la_dlarfgp(q - i,x21(i,i + 1),x21(i,i + 2),ldx21,tauq1(i))
                 s = x21(i,i + 1)
                 x21(i,i + 1) = one
                 call la_dlarf('R',p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
                 call la_dlarf('R',m - p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
                 c = sqrt(la_dnrm2(p - i,x11(i + 1,i + 1),1)**2 + la_dnrm2(m - p - i,x21(i + 1, &
                           i + 1),1)**2)
                 phi(i) = atan2(s,c)
                 call la_dorbdb5(p - i,m - p - i,q - i - 1,x11(i + 1,i + 1),1,x21(i + 1,i + 1),1,x11(i + 1, &
                           i + 2),ldx11,x21(i + 1,i + 2),ldx21,work(iorbdb5),lorbdb5,childinfo)
              end if
           end do
           return
     end subroutine la_dorbdb1
     !> QORBDB1: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
     !> M-P, or M-Q. Routines QORBDB2, QORBDB3, and QORBDB4 handle cases in
     !> which Q is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_qorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           real(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < q .or. m - p < q) then
              info = -2
           else if (q < 0 .or. m - q < q) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 2
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QORBDB1',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., q of x11 and x21
           do i = 1,q
              call la_qlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              theta(i) = atan2(x21(i,i),x11(i,i))
              c = cos(theta(i))
              s = sin(theta(i))
              x11(i,i) = one
              x21(i,i) = one
              call la_qlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
              call la_qlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
              if (i < q) then
                 call la_qrot(q - i,x11(i,i + 1),ldx11,x21(i,i + 1),ldx21,c,s)
                 call la_qlarfgp(q - i,x21(i,i + 1),x21(i,i + 2),ldx21,tauq1(i))
                 s = x21(i,i + 1)
                 x21(i,i + 1) = one
                 call la_qlarf('R',p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
                 call la_qlarf('R',m - p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
                 c = sqrt(la_qnrm2(p - i,x11(i + 1,i + 1),1)**2 + la_qnrm2(m - p - i,x21(i + 1, &
                           i + 1),1)**2)
                 phi(i) = atan2(s,c)
                 call la_qorbdb5(p - i,m - p - i,q - i - 1,x11(i + 1,i + 1),1,x21(i + 1,i + 1),1,x11(i + 1, &
                           i + 2),ldx11,x21(i + 1,i + 2),ldx21,work(iorbdb5),lorbdb5,childinfo)
              end if
           end do
           return
     end subroutine la_qorbdb1

     !> SORBDB2: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
     !> Q, or M-Q. Routines SORBDB1, SORBDB3, and SORBDB4 handle cases in
     !> which P is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_sorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           real(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < 0 .or. p > m - p) then
              info = -2
           else if (q < 0 .or. q < p .or. m - q < p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SORBDB2',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., p of x11 and x21
           do i = 1,p
              if (i > 1) then
                 call la_srot(q - i + 1,x11(i,i),ldx11,x21(i - 1,i),ldx21,c,s)
              end if
              call la_slarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              c = x11(i,i)
              x11(i,i) = one
              call la_slarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_slarf('R',m - p - i + 1,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(i,i),ldx21, &
                        work(ilarf))
              s = sqrt(la_snrm2(p - i,x11(i + 1,i),1)**2 + la_snrm2(m - p - i + 1,x21(i,i),1 &
                        )**2)
              theta(i) = atan2(s,c)
              call la_sorbdb5(p - i,m - p - i + 1,q - i,x11(i + 1,i),1,x21(i,i),1,x11(i + 1,i + 1), &
                        ldx11,x21(i,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_sscal(p - i,negone,x11(i + 1,i),1)
              call la_slarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              if (i < p) then
                 call la_slarfgp(p - i,x11(i + 1,i),x11(i + 2,i),1,taup1(i))
                 phi(i) = atan2(x11(i + 1,i),x21(i,i))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x11(i + 1,i) = one
                 call la_slarf('L',p - i,q - i,x11(i + 1,i),1,taup1(i),x11(i + 1,i + 1),ldx11, &
                           work(ilarf))
              end if
              x21(i,i) = one
              call la_slarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
           end do
           ! reduce the bottom-right portion of x21 to the identity matrix
           do i = p + 1,q
              call la_slarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              x21(i,i) = one
              call la_slarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
           end do
           return
     end subroutine la_sorbdb2
     !> DORBDB2: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
     !> Q, or M-Q. Routines DORBDB1, DORBDB3, and DORBDB4 handle cases in
     !> which P is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_dorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           real(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < 0 .or. p > m - p) then
              info = -2
           else if (q < 0 .or. q < p .or. m - q < p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DORBDB2',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., p of x11 and x21
           do i = 1,p
              if (i > 1) then
                 call la_drot(q - i + 1,x11(i,i),ldx11,x21(i - 1,i),ldx21,c,s)
              end if
              call la_dlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              c = x11(i,i)
              x11(i,i) = one
              call la_dlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_dlarf('R',m - p - i + 1,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(i,i),ldx21, &
                        work(ilarf))
              s = sqrt(la_dnrm2(p - i,x11(i + 1,i),1)**2 + la_dnrm2(m - p - i + 1,x21(i,i),1 &
                        )**2)
              theta(i) = atan2(s,c)
              call la_dorbdb5(p - i,m - p - i + 1,q - i,x11(i + 1,i),1,x21(i,i),1,x11(i + 1,i + 1), &
                        ldx11,x21(i,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_dscal(p - i,negone,x11(i + 1,i),1)
              call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              if (i < p) then
                 call la_dlarfgp(p - i,x11(i + 1,i),x11(i + 2,i),1,taup1(i))
                 phi(i) = atan2(x11(i + 1,i),x21(i,i))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x11(i + 1,i) = one
                 call la_dlarf('L',p - i,q - i,x11(i + 1,i),1,taup1(i),x11(i + 1,i + 1),ldx11, &
                           work(ilarf))
              end if
              x21(i,i) = one
              call la_dlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
           end do
           ! reduce the bottom-right portion of x21 to the identity matrix
           do i = p + 1,q
              call la_dlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              x21(i,i) = one
              call la_dlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
           end do
           return
     end subroutine la_dorbdb2
     !> QORBDB2: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
     !> Q, or M-Q. Routines QORBDB1, QORBDB3, and QORBDB4 handle cases in
     !> which P is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_qorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           real(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < 0 .or. p > m - p) then
              info = -2
           else if (q < 0 .or. q < p .or. m - q < p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QORBDB2',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., p of x11 and x21
           do i = 1,p
              if (i > 1) then
                 call la_qrot(q - i + 1,x11(i,i),ldx11,x21(i - 1,i),ldx21,c,s)
              end if
              call la_qlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              c = x11(i,i)
              x11(i,i) = one
              call la_qlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_qlarf('R',m - p - i + 1,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(i,i),ldx21, &
                        work(ilarf))
              s = sqrt(la_qnrm2(p - i,x11(i + 1,i),1)**2 + la_qnrm2(m - p - i + 1,x21(i,i),1 &
                        )**2)
              theta(i) = atan2(s,c)
              call la_qorbdb5(p - i,m - p - i + 1,q - i,x11(i + 1,i),1,x21(i,i),1,x11(i + 1,i + 1), &
                        ldx11,x21(i,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_qscal(p - i,negone,x11(i + 1,i),1)
              call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              if (i < p) then
                 call la_qlarfgp(p - i,x11(i + 1,i),x11(i + 2,i),1,taup1(i))
                 phi(i) = atan2(x11(i + 1,i),x21(i,i))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x11(i + 1,i) = one
                 call la_qlarf('L',p - i,q - i,x11(i + 1,i),1,taup1(i),x11(i + 1,i + 1),ldx11, &
                           work(ilarf))
              end if
              x21(i,i) = one
              call la_qlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
           end do
           ! reduce the bottom-right portion of x21 to the identity matrix
           do i = p + 1,q
              call la_qlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              x21(i,i) = one
              call la_qlarf('L',m - p - i + 1,q - i,x21(i,i),1,taup2(i),x21(i,i + 1),ldx21,work( &
                        ilarf))
           end do
           return
     end subroutine la_qorbdb2

     !> SORBDB3: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
     !> Q, or M-Q. Routines SORBDB1, SORBDB2, and SORBDB4 handle cases in
     !> which M-P is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_sorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           real(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (2*p < m .or. p > m) then
              info = -2
           else if (q < m - p .or. m - q < m - p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SORBDB3',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., m-p of x11 and x21
           do i = 1,m - p
              if (i > 1) then
                 call la_srot(q - i + 1,x11(i - 1,i),ldx11,x21(i,i),ldx11,c,s)
              end if
              call la_slarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              s = x21(i,i)
              x21(i,i) = one
              call la_slarf('R',p - i + 1,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i,i),ldx11, &
                        work(ilarf))
              call la_slarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              c = sqrt(la_snrm2(p - i + 1,x11(i,i),1)**2 + la_snrm2(m - p - i,x21(i + 1,i),1 &
                        )**2)
              theta(i) = atan2(s,c)
              call la_sorbdb5(p - i + 1,m - p - i,q - i,x11(i,i),1,x21(i + 1,i),1,x11(i,i + 1), &
                        ldx11,x21(i + 1,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_slarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              if (i < m - p) then
                 call la_slarfgp(m - p - i,x21(i + 1,i),x21(i + 2,i),1,taup2(i))
                 phi(i) = atan2(x21(i + 1,i),x11(i,i))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x21(i + 1,i) = one
                 call la_slarf('L',m - p - i,q - i,x21(i + 1,i),1,taup2(i),x21(i + 1,i + 1),ldx21, &
                           work(ilarf))
              end if
              x11(i,i) = one
              call la_slarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
           end do
           ! reduce the bottom-right portion of x11 to the identity matrix
           do i = m - p + 1,q
              call la_slarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              x11(i,i) = one
              call la_slarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
           end do
           return
     end subroutine la_sorbdb3
     !> DORBDB3: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
     !> Q, or M-Q. Routines DORBDB1, DORBDB2, and DORBDB4 handle cases in
     !> which M-P is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_dorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           real(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (2*p < m .or. p > m) then
              info = -2
           else if (q < m - p .or. m - q < m - p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DORBDB3',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., m-p of x11 and x21
           do i = 1,m - p
              if (i > 1) then
                 call la_drot(q - i + 1,x11(i - 1,i),ldx11,x21(i,i),ldx11,c,s)
              end if
              call la_dlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              s = x21(i,i)
              x21(i,i) = one
              call la_dlarf('R',p - i + 1,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i,i),ldx11, &
                        work(ilarf))
              call la_dlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              c = sqrt(la_dnrm2(p - i + 1,x11(i,i),1)**2 + la_dnrm2(m - p - i,x21(i + 1,i),1 &
                        )**2)
              theta(i) = atan2(s,c)
              call la_dorbdb5(p - i + 1,m - p - i,q - i,x11(i,i),1,x21(i + 1,i),1,x11(i,i + 1), &
                        ldx11,x21(i + 1,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_dlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              if (i < m - p) then
                 call la_dlarfgp(m - p - i,x21(i + 1,i),x21(i + 2,i),1,taup2(i))
                 phi(i) = atan2(x21(i + 1,i),x11(i,i))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x21(i + 1,i) = one
                 call la_dlarf('L',m - p - i,q - i,x21(i + 1,i),1,taup2(i),x21(i + 1,i + 1),ldx21, &
                           work(ilarf))
              end if
              x11(i,i) = one
              call la_dlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
           end do
           ! reduce the bottom-right portion of x11 to the identity matrix
           do i = m - p + 1,q
              call la_dlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              x11(i,i) = one
              call la_dlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
           end do
           return
     end subroutine la_dorbdb3
     !> QORBDB3: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
     !> Q, or M-Q. Routines QORBDB1, QORBDB2, and QORBDB4 handle cases in
     !> which M-P is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_qorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           real(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           real(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (2*p < m .or. p > m) then
              info = -2
           else if (q < m - p .or. m - q < m - p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QORBDB3',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., m-p of x11 and x21
           do i = 1,m - p
              if (i > 1) then
                 call la_qrot(q - i + 1,x11(i - 1,i),ldx11,x21(i,i),ldx11,c,s)
              end if
              call la_qlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              s = x21(i,i)
              x21(i,i) = one
              call la_qlarf('R',p - i + 1,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i,i),ldx11, &
                        work(ilarf))
              call la_qlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              c = sqrt(la_qnrm2(p - i + 1,x11(i,i),1)**2 + la_qnrm2(m - p - i,x21(i + 1,i),1 &
                        )**2)
              theta(i) = atan2(s,c)
              call la_qorbdb5(p - i + 1,m - p - i,q - i,x11(i,i),1,x21(i + 1,i),1,x11(i,i + 1), &
                        ldx11,x21(i + 1,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_qlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              if (i < m - p) then
                 call la_qlarfgp(m - p - i,x21(i + 1,i),x21(i + 2,i),1,taup2(i))
                 phi(i) = atan2(x21(i + 1,i),x11(i,i))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x21(i + 1,i) = one
                 call la_qlarf('L',m - p - i,q - i,x21(i + 1,i),1,taup2(i),x21(i + 1,i + 1),ldx21, &
                           work(ilarf))
              end if
              x11(i,i) = one
              call la_qlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
           end do
           ! reduce the bottom-right portion of x11 to the identity matrix
           do i = m - p + 1,q
              call la_qlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              x11(i,i) = one
              call la_qlarf('L',p - i + 1,q - i,x11(i,i),1,taup1(i),x11(i,i + 1),ldx11,work( &
                        ilarf))
           end do
           return
     end subroutine la_qorbdb3

     !> SORBDB4: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
     !> M-P, or Q. Routines SORBDB1, SORBDB2, and SORBDB3 handle cases in
     !> which M-Q is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_sorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               phantom,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           real(sp),intent(out) :: phantom(*),taup1(*),taup2(*),tauq1(*),work(*)
           real(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,j,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < m - q .or. m - p < m - q) then
              info = -2
           else if (q < m - q .or. q > m) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(q - 1,p - 1,m - p - 1)
              iorbdb5 = 2
              lorbdb5 = q
              lworkopt = ilarf + llarf - 1
              lworkopt = max(lworkopt,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SORBDB4',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., m-q of x11 and x21
           do i = 1,m - q
              if (i == 1) then
                 do j = 1,m
                    phantom(j) = zero
                 end do
                 call la_sorbdb5(p,m - p,q,phantom(1),1,phantom(p + 1),1,x11,ldx11,x21, &
                           ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_sscal(p,negone,phantom(1),1)
                 call la_slarfgp(p,phantom(1),phantom(2),1,taup1(1))
                 call la_slarfgp(m - p,phantom(p + 1),phantom(p + 2),1,taup2(1))
                 theta(i) = atan2(phantom(1),phantom(p + 1))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 phantom(1) = one
                 phantom(p + 1) = one
                 call la_slarf('L',p,q,phantom(1),1,taup1(1),x11,ldx11,work(ilarf))

                 call la_slarf('L',m - p,q,phantom(p + 1),1,taup2(1),x21,ldx21,work(ilarf) &
                            )
              else
                 call la_sorbdb5(p - i + 1,m - p - i + 1,q - i + 1,x11(i,i - 1),1,x21(i,i - 1),1,x11(i,i) &
                           ,ldx11,x21(i,i),ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_sscal(p - i + 1,negone,x11(i,i - 1),1)
                 call la_slarfgp(p - i + 1,x11(i,i - 1),x11(i + 1,i - 1),1,taup1(i))
                 call la_slarfgp(m - p - i + 1,x21(i,i - 1),x21(i + 1,i - 1),1,taup2(i))
                 theta(i) = atan2(x11(i,i - 1),x21(i,i - 1))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 x11(i,i - 1) = one
                 x21(i,i - 1) = one
                 call la_slarf('L',p - i + 1,q - i + 1,x11(i,i - 1),1,taup1(i),x11(i,i),ldx11, &
                           work(ilarf))
                 call la_slarf('L',m - p - i + 1,q - i + 1,x21(i,i - 1),1,taup2(i),x21(i,i),ldx21, &
                           work(ilarf))
              end if
              call la_srot(q - i + 1,x11(i,i),ldx11,x21(i,i),ldx21,s,-c)
              call la_slarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              c = x21(i,i)
              x21(i,i) = one
              call la_slarf('R',p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_slarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              if (i < m - q) then
                 s = sqrt(la_snrm2(p - i,x11(i + 1,i),1)**2 + la_snrm2(m - p - i,x21(i + 1,i), &
                            1)**2)
                 phi(i) = atan2(s,c)
              end if
           end do
           ! reduce the bottom-right portion of x11 to [ i 0 ]
           do i = m - q + 1,p
              call la_slarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              x11(i,i) = one
              call la_slarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_slarf('R',q - p,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(m - q + 1,i),ldx21, &
                        work(ilarf))
           end do
           ! reduce the bottom-right portion of x21 to [ 0 i ]
           do i = p + 1,q
              call la_slarfgp(q - i + 1,x21(m - q + i - p,i),x21(m - q + i - p,i + 1),ldx21,tauq1(i))

              x21(m - q + i - p,i) = one
              call la_slarf('R',q - i,q - i + 1,x21(m - q + i - p,i),ldx21,tauq1(i),x21(m - q + i - p + 1,i) &
                        ,ldx21,work(ilarf))
           end do
           return
     end subroutine la_sorbdb4
     !> DORBDB4: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
     !> M-P, or Q. Routines DORBDB1, DORBDB2, and DORBDB3 handle cases in
     !> which M-Q is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_dorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               phantom,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           real(dp),intent(out) :: phantom(*),taup1(*),taup2(*),tauq1(*),work(*)
           real(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,j,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < m - q .or. m - p < m - q) then
              info = -2
           else if (q < m - q .or. q > m) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(q - 1,p - 1,m - p - 1)
              iorbdb5 = 2
              lorbdb5 = q
              lworkopt = ilarf + llarf - 1
              lworkopt = max(lworkopt,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DORBDB4',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., m-q of x11 and x21
           do i = 1,m - q
              if (i == 1) then
                 do j = 1,m
                    phantom(j) = zero
                 end do
                 call la_dorbdb5(p,m - p,q,phantom(1),1,phantom(p + 1),1,x11,ldx11,x21, &
                           ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_dscal(p,negone,phantom(1),1)
                 call la_dlarfgp(p,phantom(1),phantom(2),1,taup1(1))
                 call la_dlarfgp(m - p,phantom(p + 1),phantom(p + 2),1,taup2(1))
                 theta(i) = atan2(phantom(1),phantom(p + 1))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 phantom(1) = one
                 phantom(p + 1) = one
                 call la_dlarf('L',p,q,phantom(1),1,taup1(1),x11,ldx11,work(ilarf))

                 call la_dlarf('L',m - p,q,phantom(p + 1),1,taup2(1),x21,ldx21,work(ilarf) &
                            )
              else
                 call la_dorbdb5(p - i + 1,m - p - i + 1,q - i + 1,x11(i,i - 1),1,x21(i,i - 1),1,x11(i,i) &
                           ,ldx11,x21(i,i),ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_dscal(p - i + 1,negone,x11(i,i - 1),1)
                 call la_dlarfgp(p - i + 1,x11(i,i - 1),x11(i + 1,i - 1),1,taup1(i))
                 call la_dlarfgp(m - p - i + 1,x21(i,i - 1),x21(i + 1,i - 1),1,taup2(i))
                 theta(i) = atan2(x11(i,i - 1),x21(i,i - 1))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 x11(i,i - 1) = one
                 x21(i,i - 1) = one
                 call la_dlarf('L',p - i + 1,q - i + 1,x11(i,i - 1),1,taup1(i),x11(i,i),ldx11, &
                           work(ilarf))
                 call la_dlarf('L',m - p - i + 1,q - i + 1,x21(i,i - 1),1,taup2(i),x21(i,i),ldx21, &
                           work(ilarf))
              end if
              call la_drot(q - i + 1,x11(i,i),ldx11,x21(i,i),ldx21,s,-c)
              call la_dlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              c = x21(i,i)
              x21(i,i) = one
              call la_dlarf('R',p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_dlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              if (i < m - q) then
                 s = sqrt(la_dnrm2(p - i,x11(i + 1,i),1)**2 + la_dnrm2(m - p - i,x21(i + 1,i), &
                            1)**2)
                 phi(i) = atan2(s,c)
              end if
           end do
           ! reduce the bottom-right portion of x11 to [ i 0 ]
           do i = m - q + 1,p
              call la_dlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              x11(i,i) = one
              call la_dlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_dlarf('R',q - p,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(m - q + 1,i),ldx21, &
                        work(ilarf))
           end do
           ! reduce the bottom-right portion of x21 to [ 0 i ]
           do i = p + 1,q
              call la_dlarfgp(q - i + 1,x21(m - q + i - p,i),x21(m - q + i - p,i + 1),ldx21,tauq1(i))

              x21(m - q + i - p,i) = one
              call la_dlarf('R',q - i,q - i + 1,x21(m - q + i - p,i),ldx21,tauq1(i),x21(m - q + i - p + 1,i) &
                        ,ldx21,work(ilarf))
           end do
           return
     end subroutine la_dorbdb4
     !> QORBDB4: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
     !> M-P, or Q. Routines QORBDB1, QORBDB2, and QORBDB3 handle cases in
     !> which M-Q is not the minimum dimension.
     !> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_qorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               phantom,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           real(qp),intent(out) :: phantom(*),taup1(*),taup2(*),tauq1(*),work(*)
           real(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,j,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < m - q .or. m - p < m - q) then
              info = -2
           else if (q < m - q .or. q > m) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(q - 1,p - 1,m - p - 1)
              iorbdb5 = 2
              lorbdb5 = q
              lworkopt = ilarf + llarf - 1
              lworkopt = max(lworkopt,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QORBDB4',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., m-q of x11 and x21
           do i = 1,m - q
              if (i == 1) then
                 do j = 1,m
                    phantom(j) = zero
                 end do
                 call la_qorbdb5(p,m - p,q,phantom(1),1,phantom(p + 1),1,x11,ldx11,x21, &
                           ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_qscal(p,negone,phantom(1),1)
                 call la_qlarfgp(p,phantom(1),phantom(2),1,taup1(1))
                 call la_qlarfgp(m - p,phantom(p + 1),phantom(p + 2),1,taup2(1))
                 theta(i) = atan2(phantom(1),phantom(p + 1))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 phantom(1) = one
                 phantom(p + 1) = one
                 call la_qlarf('L',p,q,phantom(1),1,taup1(1),x11,ldx11,work(ilarf))

                 call la_qlarf('L',m - p,q,phantom(p + 1),1,taup2(1),x21,ldx21,work(ilarf) &
                            )
              else
                 call la_qorbdb5(p - i + 1,m - p - i + 1,q - i + 1,x11(i,i - 1),1,x21(i,i - 1),1,x11(i,i) &
                           ,ldx11,x21(i,i),ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_qscal(p - i + 1,negone,x11(i,i - 1),1)
                 call la_qlarfgp(p - i + 1,x11(i,i - 1),x11(i + 1,i - 1),1,taup1(i))
                 call la_qlarfgp(m - p - i + 1,x21(i,i - 1),x21(i + 1,i - 1),1,taup2(i))
                 theta(i) = atan2(x11(i,i - 1),x21(i,i - 1))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 x11(i,i - 1) = one
                 x21(i,i - 1) = one
                 call la_qlarf('L',p - i + 1,q - i + 1,x11(i,i - 1),1,taup1(i),x11(i,i),ldx11, &
                           work(ilarf))
                 call la_qlarf('L',m - p - i + 1,q - i + 1,x21(i,i - 1),1,taup2(i),x21(i,i),ldx21, &
                           work(ilarf))
              end if
              call la_qrot(q - i + 1,x11(i,i),ldx11,x21(i,i),ldx21,s,-c)
              call la_qlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              c = x21(i,i)
              x21(i,i) = one
              call la_qlarf('R',p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_qlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              if (i < m - q) then
                 s = sqrt(la_qnrm2(p - i,x11(i + 1,i),1)**2 + la_qnrm2(m - p - i,x21(i + 1,i), &
                            1)**2)
                 phi(i) = atan2(s,c)
              end if
           end do
           ! reduce the bottom-right portion of x11 to [ i 0 ]
           do i = m - q + 1,p
              call la_qlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              x11(i,i) = one
              call la_qlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_qlarf('R',q - p,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(m - q + 1,i),ldx21, &
                        work(ilarf))
           end do
           ! reduce the bottom-right portion of x21 to [ 0 i ]
           do i = p + 1,q
              call la_qlarfgp(q - i + 1,x21(m - q + i - p,i),x21(m - q + i - p,i + 1),ldx21,tauq1(i))

              x21(m - q + i - p,i) = one
              call la_qlarf('R',q - i,q - i + 1,x21(m - q + i - p,i),ldx21,tauq1(i),x21(m - q + i - p + 1,i) &
                        ,ldx21,work(ilarf))
           end do
           return
     end subroutine la_qorbdb4

     !> SORCSD2BY1: computes the CS decomposition of an M-by-Q matrix X with
     !> orthonormal columns that has been partitioned into a 2-by-1 block
     !> structure:
     !> [  I1 0  0 ]
     !> [  0  C  0 ]
     !> [ X11 ]   [ U1 |    ] [  0  0  0 ]
     !> X = [-----] = [---------] [----------] V1**T .
     !> [ X21 ]   [    | U2 ] [  0  0  0 ]
     !> [  0  S  0 ]
     !> [  0  0  I2]
     !> X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P,
     !> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     !> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     !> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     !> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).

     subroutine la_sorcsd2by1(jobu1,jobu2,jobv1t,m,p,q,x11,ldx11,x21,ldx21,theta, &
               u1,ldu1,u2,ldu2,v1t,ldv1t,work,lwork,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,lwork,ldx11,ldx21,m,p,q
           ! Array Arguments
           real(sp),intent(out) :: theta(*)
           real(sp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),work(*)
           real(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
           integer(ilp),intent(out) :: iwork(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,j,lbbcsd,lorbdb, &
           lorglq,lorglqmin,lorglqopt,lorgqr,lorgqrmin,lorgqropt,lworkmin,lworkopt, &
                     r
           logical(lk) :: lquery,wantu1,wantu2,wantv1t
           ! Local Arrays
           real(sp) :: dum1(1),dum2(1,1)
           ! Intrinsic Function
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           lquery = lwork == -1
           if (m < 0) then
              info = -4
           else if (p < 0 .or. p > m) then
              info = -5
           else if (q < 0 .or. q > m) then
              info = -6
           else if (ldx11 < max(1,p)) then
              info = -8
           else if (ldx21 < max(1,m - p)) then
              info = -10
           else if (wantu1 .and. ldu1 < max(1,p)) then
              info = -13
           else if (wantu2 .and. ldu2 < max(1,m - p)) then
              info = -15
           else if (wantv1t .and. ldv1t < max(1,q)) then
              info = -17
           end if
           r = min(p,m - p,q,m - q)
           ! compute workspace
             ! work layout:
           ! |-------------------------------------------------------|
           ! | lworkopt (1)                                          |
           ! |-------------------------------------------------------|
           ! | phi (max(1,r-1))                                      |
           ! |-------------------------------------------------------|
           ! | taup1 (max(1,p))                        | b11d (r)    |
           ! | taup2 (max(1,m-p))                      | b11e (r-1)  |
           ! | tauq1 (max(1,q))                        | b12d (r)    |
           ! |-----------------------------------------| b12e (r-1)  |
           ! | la_sorbdb work | la_sorgqr work | la_sorglq work | b21d (r)    |
           ! |             |             |             | b21e (r-1)  |
           ! |             |             |             | b22d (r)    |
           ! |             |             |             | b22e (r-1)  |
           ! |             |             |             | la_sbbcsd work |
           ! |-------------------------------------------------------|
           if (info == 0) then
              iphi = 2
              ib11d = iphi + max(1,r - 1)
              ib11e = ib11d + max(1,r)
              ib12d = ib11e + max(1,r - 1)
              ib12e = ib12d + max(1,r)
              ib21d = ib12e + max(1,r - 1)
              ib21e = ib21d + max(1,r)
              ib22d = ib21e + max(1,r - 1)
              ib22e = ib22d + max(1,r)
              ibbcsd = ib22e + max(1,r - 1)
              itaup1 = iphi + max(1,r - 1)
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              iorbdb = itauq1 + max(1,q)
              iorgqr = itauq1 + max(1,q)
              iorglq = itauq1 + max(1,q)
              lorgqrmin = 1
              lorgqropt = 1
              lorglqmin = 1
              lorglqopt = 1
              if (r == q) then
                 call la_sorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work,-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_sorgqr(p,p,q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_sorgqr(m - p,m - p,q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_sorglq(q - 1,q - 1,q - 1,v1t,ldv1t,dum1,work(1),-1,childinfo)

                    lorglqmin = max(lorglqmin,q - 1)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_sbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,dum1,u1, &
                 ldu1,u2,ldu2,v1t,ldv1t,dum2,1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else if (r == p) then
                 call la_sorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_sorgqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,dum1,work(1),-1,childinfo &
                              )
                    lorgqrmin = max(lorgqrmin,p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_sorgqr(m - p,m - p,q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_sorglq(q,q,r,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_sbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,dum1,v1t, &
                 ldv1t,dum2,1,u1,ldu1,u2,ldu2,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else if (r == m - p) then
                 call la_sorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_sorgqr(p,p,q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_sorgqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,dum1,work(1),-1, &
                              childinfo)
                    lorgqrmin = max(lorgqrmin,m - p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_sorglq(q,q,r,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_sbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,dum1, &
                 dum2,1,v1t,ldv1t,u2,ldu2,u1,ldu1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else
                 call la_sorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,dum1,work(1),-1,childinfo)
                 lorbdb = m + int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_sorgqr(p,p,m - q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_sorgqr(m - p,m - p,m - q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_sorglq(q,q,q,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_sbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,dum1,u2, &
                 ldu2,u1,ldu1,dum2,1,v1t,ldv1t,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              end if
              lworkmin = max(iorbdb + lorbdb - 1,iorgqr + lorgqrmin - 1,iorglq + lorglqmin - 1,ibbcsd + lbbcsd - &
                        1)
              lworkopt = max(iorbdb + lorbdb - 1,iorgqr + lorgqropt - 1,iorglq + lorglqopt - 1,ibbcsd + lbbcsd - &
                        1)
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -19
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SORCSD2BY1',-info)
              return
           else if (lquery) then
              return
           end if
           lorgqr = lwork - iorgqr + 1
           lorglq = lwork - iorglq + 1
           ! handle four cases separately: r = q, r = p, r = m-p, and r = m-q,
           ! in which r = min(p,m-p,q,m-q)
           if (r == q) then
              ! case 1: r = q
              ! simultaneously bidiagonalize x11 and x21
              call la_sorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_slacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_sorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_slacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_sorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_slacpy('U',q - 1,q - 1,x21(1,2),ldx21,v1t(2,2),ldv1t)
                 call la_sorglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglq,childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_sbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,work(iphi),u1, &
              ldu1,u2,ldu2,v1t,ldv1t,dum2,1,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place zero submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_slapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == p) then
              ! case 2: r = p
              ! simultaneously bidiagonalize x11 and x21
              call la_sorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 u1(1,1) = one
                 do j = 2,p
                    u1(1,j) = zero
                    u1(j,1) = zero
                 end do
                 call la_slacpy('L',p - 1,p - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_sorgqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,work(itaup1),work(iorgqr), &
                           lorgqr,childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_slacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_sorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_slacpy('U',p,q,x11,ldx11,v1t,ldv1t)
                 call la_sorglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_sbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,work(iphi),v1t, &
              ldv1t,dum1,1,u1,ldu1,u2,ldu2,work(ib11d),work(ib11e),work(ib12d),work(ib12e) &
              ,work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd,childinfo &
                        )
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_slapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == m - p) then
              ! case 3: r = m-p
              ! simultaneously bidiagonalize x11 and x21
              call la_sorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_slacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_sorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 u2(1,1) = one
                 do j = 2,m - p
                    u2(1,j) = zero
                    u2(j,1) = zero
                 end do
                 call la_slacpy('L',m - p - 1,m - p - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_sorgqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,work(itaup2),work(iorgqr) &
                           ,lorgqr,childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_slacpy('U',m - p,q,x21,ldx21,v1t,ldv1t)
                 call la_sorglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_sbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,work(iphi), &
              dum1,1,v1t,ldv1t,u2,ldu2,u1,ldu1,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > r) then
                 do i = 1,r
                    iwork(i) = q - r + i
                 end do
                 do i = r + 1,q
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_slapmt(.false.,p,q,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_slapmr(.false.,q,q,v1t,ldv1t,iwork)
                 end if
              end if
           else
              ! case 4: r = m-q
              ! simultaneously bidiagonalize x11 and x21
              call la_sorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
              ,work(itaup2),work(itauq1),work(iorbdb),work(iorbdb + m),lorbdb - m,childinfo)

              ! accumulate householder reflectors
              if (wantu2 .and. m - p > 0) then
                 call la_scopy(m - p,work(iorbdb + p),1,u2,1)
              end if
              if (wantu1 .and. p > 0) then
                 call la_scopy(p,work(iorbdb),1,u1,1)
                 do j = 2,p
                    u1(1,j) = zero
                 end do
                 call la_slacpy('L',p - 1,m - q - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_sorgqr(p,p,m - q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 do j = 2,m - p
                    u2(1,j) = zero
                 end do
                 call la_slacpy('L',m - p - 1,m - q - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_sorgqr(m - p,m - p,m - q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_slacpy('U',m - q,q,x21,ldx21,v1t,ldv1t)
                 call la_slacpy('U',p - (m - q),q - (m - q),x11(m - q + 1,m - q + 1),ldx11,v1t(m - q + 1,m - q + &
                           1),ldv1t)
                 call la_slacpy('U',-p + q,q - p,x21(m - q + 1,p + 1),ldx21,v1t(p + 1,p + 1),ldv1t)

                 call la_sorglq(q,q,q,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_sbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,work(iphi), &
              u2,ldu2,u1,ldu1,dum1,1,v1t,ldv1t,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (p > r) then
                 do i = 1,r
                    iwork(i) = p - r + i
                 end do
                 do i = r + 1,p
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_slapmt(.false.,p,p,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_slapmr(.false.,p,q,v1t,ldv1t,iwork)
                 end if
              end if
           end if
           return
     end subroutine la_sorcsd2by1
     !> DORCSD2BY1: computes the CS decomposition of an M-by-Q matrix X with
     !> orthonormal columns that has been partitioned into a 2-by-1 block
     !> structure:
     !> [  I1 0  0 ]
     !> [  0  C  0 ]
     !> [ X11 ]   [ U1 |    ] [  0  0  0 ]
     !> X = [-----] = [---------] [----------] V1**T .
     !> [ X21 ]   [    | U2 ] [  0  0  0 ]
     !> [  0  S  0 ]
     !> [  0  0  I2]
     !> X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P,
     !> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     !> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     !> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     !> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).

     subroutine la_dorcsd2by1(jobu1,jobu2,jobv1t,m,p,q,x11,ldx11,x21,ldx21,theta, &
               u1,ldu1,u2,ldu2,v1t,ldv1t,work,lwork,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine (3.5.0_dp) --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,lwork,ldx11,ldx21,m,p,q
           ! Array Arguments
           real(dp),intent(out) :: theta(*)
           real(dp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),work(*)
           real(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
           integer(ilp),intent(out) :: iwork(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,j,lbbcsd,lorbdb, &
           lorglq,lorglqmin,lorglqopt,lorgqr,lorgqrmin,lorgqropt,lworkmin,lworkopt, &
                     r
           logical(lk) :: lquery,wantu1,wantu2,wantv1t
           ! Local Arrays
           real(dp) :: dum1(1),dum2(1,1)
           ! Intrinsic Function
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           lquery = lwork == -1
           if (m < 0) then
              info = -4
           else if (p < 0 .or. p > m) then
              info = -5
           else if (q < 0 .or. q > m) then
              info = -6
           else if (ldx11 < max(1,p)) then
              info = -8
           else if (ldx21 < max(1,m - p)) then
              info = -10
           else if (wantu1 .and. ldu1 < max(1,p)) then
              info = -13
           else if (wantu2 .and. ldu2 < max(1,m - p)) then
              info = -15
           else if (wantv1t .and. ldv1t < max(1,q)) then
              info = -17
           end if
           r = min(p,m - p,q,m - q)
           ! compute workspace
             ! work layout:
           ! |-------------------------------------------------------|
           ! | lworkopt (1)                                          |
           ! |-------------------------------------------------------|
           ! | phi (max(1,r-1))                                      |
           ! |-------------------------------------------------------|
           ! | taup1 (max(1,p))                        | b11d (r)    |
           ! | taup2 (max(1,m-p))                      | b11e (r-1)  |
           ! | tauq1 (max(1,q))                        | b12d (r)    |
           ! |-----------------------------------------| b12e (r-1)  |
           ! | la_dorbdb work | la_dorgqr work | la_dorglq work | b21d (r)    |
           ! |             |             |             | b21e (r-1)  |
           ! |             |             |             | b22d (r)    |
           ! |             |             |             | b22e (r-1)  |
           ! |             |             |             | la_dbbcsd work |
           ! |-------------------------------------------------------|
           if (info == 0) then
              iphi = 2
              ib11d = iphi + max(1,r - 1)
              ib11e = ib11d + max(1,r)
              ib12d = ib11e + max(1,r - 1)
              ib12e = ib12d + max(1,r)
              ib21d = ib12e + max(1,r - 1)
              ib21e = ib21d + max(1,r)
              ib22d = ib21e + max(1,r - 1)
              ib22e = ib22d + max(1,r)
              ibbcsd = ib22e + max(1,r - 1)
              itaup1 = iphi + max(1,r - 1)
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              iorbdb = itauq1 + max(1,q)
              iorgqr = itauq1 + max(1,q)
              iorglq = itauq1 + max(1,q)
              lorgqrmin = 1
              lorgqropt = 1
              lorglqmin = 1
              lorglqopt = 1
              if (r == q) then
                 call la_dorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work,-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_dorgqr(p,p,q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_dorgqr(m - p,m - p,q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_dorglq(q - 1,q - 1,q - 1,v1t,ldv1t,dum1,work(1),-1,childinfo)

                    lorglqmin = max(lorglqmin,q - 1)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_dbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,dum1,u1, &
                 ldu1,u2,ldu2,v1t,ldv1t,dum2,1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                            work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else if (r == p) then
                 call la_dorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_dorgqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,dum1,work(1),-1,childinfo &
                              )
                    lorgqrmin = max(lorgqrmin,p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_dorgqr(m - p,m - p,q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_dorglq(q,q,r,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_dbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,dum1,v1t, &
                 ldv1t,dum2,1,u1,ldu1,u2,ldu2,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else if (r == m - p) then
                 call la_dorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_dorgqr(p,p,q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_dorgqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,dum1,work(1),-1, &
                              childinfo)
                    lorgqrmin = max(lorgqrmin,m - p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_dorglq(q,q,r,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_dbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,dum1, &
                 dum2,1,v1t,ldv1t,u2,ldu2,u1,ldu1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else
                 call la_dorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,dum1,work(1),-1,childinfo)
                 lorbdb = m + int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_dorgqr(p,p,m - q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_dorgqr(m - p,m - p,m - q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_dorglq(q,q,q,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_dbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,dum1,u2, &
                 ldu2,u1,ldu1,dum2,1,v1t,ldv1t,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                            work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              end if
              lworkmin = max(iorbdb + lorbdb - 1,iorgqr + lorgqrmin - 1,iorglq + lorglqmin - 1,ibbcsd + lbbcsd - &
                        1)
              lworkopt = max(iorbdb + lorbdb - 1,iorgqr + lorgqropt - 1,iorglq + lorglqopt - 1,ibbcsd + lbbcsd - &
                        1)
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -19
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DORCSD2BY1',-info)
              return
           else if (lquery) then
              return
           end if
           lorgqr = lwork - iorgqr + 1
           lorglq = lwork - iorglq + 1
           ! handle four cases separately: r = q, r = p, r = m-p, and r = m-q,
           ! in which r = min(p,m-p,q,m-q)
           if (r == q) then
              ! case 1: r = q
              ! simultaneously bidiagonalize x11 and x21
              call la_dorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_dlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_dorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_dlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_dorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_dlacpy('U',q - 1,q - 1,x21(1,2),ldx21,v1t(2,2),ldv1t)
                 call la_dorglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglq,childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_dbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,work(iphi),u1, &
              ldu1,u2,ldu2,v1t,ldv1t,dum2,1,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place zero submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_dlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == p) then
              ! case 2: r = p
              ! simultaneously bidiagonalize x11 and x21
              call la_dorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 u1(1,1) = one
                 do j = 2,p
                    u1(1,j) = zero
                    u1(j,1) = zero
                 end do
                 call la_dlacpy('L',p - 1,p - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_dorgqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,work(itaup1),work(iorgqr), &
                           lorgqr,childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_dlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_dorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_dlacpy('U',p,q,x11,ldx11,v1t,ldv1t)
                 call la_dorglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_dbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,work(iphi),v1t, &
              ldv1t,dum2,1,u1,ldu1,u2,ldu2,work(ib11d),work(ib11e),work(ib12d),work(ib12e) &
              ,work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd,childinfo &
                        )
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_dlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == m - p) then
              ! case 3: r = m-p
              ! simultaneously bidiagonalize x11 and x21
              call la_dorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_dlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_dorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 u2(1,1) = one
                 do j = 2,m - p
                    u2(1,j) = zero
                    u2(j,1) = zero
                 end do
                 call la_dlacpy('L',m - p - 1,m - p - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_dorgqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,work(itaup2),work(iorgqr) &
                           ,lorgqr,childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_dlacpy('U',m - p,q,x21,ldx21,v1t,ldv1t)
                 call la_dorglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_dbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,work(iphi), &
              dum2,1,v1t,ldv1t,u2,ldu2,u1,ldu1,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > r) then
                 do i = 1,r
                    iwork(i) = q - r + i
                 end do
                 do i = r + 1,q
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_dlapmt(.false.,p,q,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_dlapmr(.false.,q,q,v1t,ldv1t,iwork)
                 end if
              end if
           else
              ! case 4: r = m-q
              ! simultaneously bidiagonalize x11 and x21
              call la_dorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
              ,work(itaup2),work(itauq1),work(iorbdb),work(iorbdb + m),lorbdb - m,childinfo)

              ! accumulate householder reflectors
              if (wantu2 .and. m - p > 0) then
                 call la_dcopy(m - p,work(iorbdb + p),1,u2,1)
              end if
              if (wantu1 .and. p > 0) then
                 call la_dcopy(p,work(iorbdb),1,u1,1)
                 do j = 2,p
                    u1(1,j) = zero
                 end do
                 call la_dlacpy('L',p - 1,m - q - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_dorgqr(p,p,m - q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 do j = 2,m - p
                    u2(1,j) = zero
                 end do
                 call la_dlacpy('L',m - p - 1,m - q - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_dorgqr(m - p,m - p,m - q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_dlacpy('U',m - q,q,x21,ldx21,v1t,ldv1t)
                 call la_dlacpy('U',p - (m - q),q - (m - q),x11(m - q + 1,m - q + 1),ldx11,v1t(m - q + 1,m - q + &
                           1),ldv1t)
                 call la_dlacpy('U',-p + q,q - p,x21(m - q + 1,p + 1),ldx21,v1t(p + 1,p + 1),ldv1t)

                 call la_dorglq(q,q,q,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_dbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,work(iphi), &
              u2,ldu2,u1,ldu1,dum2,1,v1t,ldv1t,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (p > r) then
                 do i = 1,r
                    iwork(i) = p - r + i
                 end do
                 do i = r + 1,p
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_dlapmt(.false.,p,p,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_dlapmr(.false.,p,q,v1t,ldv1t,iwork)
                 end if
              end if
           end if
           return
     end subroutine la_dorcsd2by1
     !> QORCSD2BY1: computes the CS decomposition of an M-by-Q matrix X with
     !> orthonormal columns that has been partitioned into a 2-by-1 block
     !> structure:
     !> [  I1 0  0 ]
     !> [  0  C  0 ]
     !> [ X11 ]   [ U1 |    ] [  0  0  0 ]
     !> X = [-----] = [---------] [----------] V1**T .
     !> [ X21 ]   [    | U2 ] [  0  0  0 ]
     !> [  0  S  0 ]
     !> [  0  0  I2]
     !> X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P,
     !> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     !> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     !> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     !> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).

     subroutine la_qorcsd2by1(jobu1,jobu2,jobv1t,m,p,q,x11,ldx11,x21,ldx21,theta, &
               u1,ldu1,u2,ldu2,v1t,ldv1t,work,lwork,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine (3.5.0_qp) --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,lwork,ldx11,ldx21,m,p,q
           ! Array Arguments
           real(qp),intent(out) :: theta(*)
           real(qp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),work(*)
           real(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
           integer(ilp),intent(out) :: iwork(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,j,lbbcsd,lorbdb, &
           lorglq,lorglqmin,lorglqopt,lorgqr,lorgqrmin,lorgqropt,lworkmin,lworkopt, &
                     r
           logical(lk) :: lquery,wantu1,wantu2,wantv1t
           ! Local Arrays
           real(qp) :: dum1(1),dum2(1,1)
           ! Intrinsic Function
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           lquery = lwork == -1
           if (m < 0) then
              info = -4
           else if (p < 0 .or. p > m) then
              info = -5
           else if (q < 0 .or. q > m) then
              info = -6
           else if (ldx11 < max(1,p)) then
              info = -8
           else if (ldx21 < max(1,m - p)) then
              info = -10
           else if (wantu1 .and. ldu1 < max(1,p)) then
              info = -13
           else if (wantu2 .and. ldu2 < max(1,m - p)) then
              info = -15
           else if (wantv1t .and. ldv1t < max(1,q)) then
              info = -17
           end if
           r = min(p,m - p,q,m - q)
           ! compute workspace
             ! work layout:
           ! |-------------------------------------------------------|
           ! | lworkopt (1)                                          |
           ! |-------------------------------------------------------|
           ! | phi (max(1,r-1))                                      |
           ! |-------------------------------------------------------|
           ! | taup1 (max(1,p))                        | b11d (r)    |
           ! | taup2 (max(1,m-p))                      | b11e (r-1)  |
           ! | tauq1 (max(1,q))                        | b12d (r)    |
           ! |-----------------------------------------| b12e (r-1)  |
           ! | la_qorbdb work | la_qorgqr work | la_qorglq work | b21d (r)    |
           ! |             |             |             | b21e (r-1)  |
           ! |             |             |             | b22d (r)    |
           ! |             |             |             | b22e (r-1)  |
           ! |             |             |             | la_qbbcsd work |
           ! |-------------------------------------------------------|
           if (info == 0) then
              iphi = 2
              ib11d = iphi + max(1,r - 1)
              ib11e = ib11d + max(1,r)
              ib12d = ib11e + max(1,r - 1)
              ib12e = ib12d + max(1,r)
              ib21d = ib12e + max(1,r - 1)
              ib21e = ib21d + max(1,r)
              ib22d = ib21e + max(1,r - 1)
              ib22e = ib22d + max(1,r)
              ibbcsd = ib22e + max(1,r - 1)
              itaup1 = iphi + max(1,r - 1)
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              iorbdb = itauq1 + max(1,q)
              iorgqr = itauq1 + max(1,q)
              iorglq = itauq1 + max(1,q)
              lorgqrmin = 1
              lorgqropt = 1
              lorglqmin = 1
              lorglqopt = 1
              if (r == q) then
                 call la_qorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work,-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_qorgqr(p,p,q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_qorgqr(m - p,m - p,q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_qorglq(q - 1,q - 1,q - 1,v1t,ldv1t,dum1,work(1),-1,childinfo)

                    lorglqmin = max(lorglqmin,q - 1)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_qbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,dum1,u1, &
                 ldu1,u2,ldu2,v1t,ldv1t,dum2,1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                            work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else if (r == p) then
                 call la_qorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_qorgqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,dum1,work(1),-1,childinfo &
                              )
                    lorgqrmin = max(lorgqrmin,p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_qorgqr(m - p,m - p,q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_qorglq(q,q,r,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_qbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,dum1,v1t, &
                 ldv1t,dum2,1,u1,ldu1,u2,ldu2,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else if (r == m - p) then
                 call la_qorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_qorgqr(p,p,q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_qorgqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,dum1,work(1),-1, &
                              childinfo)
                    lorgqrmin = max(lorgqrmin,m - p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_qorglq(q,q,r,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_qbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,dum1, &
                 dum2,1,v1t,ldv1t,u2,ldu2,u1,ldu1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                           dum1,work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              else
                 call la_qorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,dum1,dum1,dum1, &
                           dum1,dum1,work(1),-1,childinfo)
                 lorbdb = m + int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_qorgqr(p,p,m - q,u1,ldu1,dum1,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_qorgqr(m - p,m - p,m - q,u2,ldu2,dum1,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_qorglq(q,q,q,v1t,ldv1t,dum1,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_qbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,dum1,u2, &
                 ldu2,u1,ldu1,dum2,1,v1t,ldv1t,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1, &
                            work(1),-1,childinfo)
                 lbbcsd = int(work(1),KIND=ilp)
              end if
              lworkmin = max(iorbdb + lorbdb - 1,iorgqr + lorgqrmin - 1,iorglq + lorglqmin - 1,ibbcsd + lbbcsd - &
                        1)
              lworkopt = max(iorbdb + lorbdb - 1,iorgqr + lorgqropt - 1,iorglq + lorglqopt - 1,ibbcsd + lbbcsd - &
                        1)
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -19
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QORCSD2BY1',-info)
              return
           else if (lquery) then
              return
           end if
           lorgqr = lwork - iorgqr + 1
           lorglq = lwork - iorglq + 1
           ! handle four cases separately: r = q, r = p, r = m-p, and r = m-q,
           ! in which r = min(p,m-p,q,m-q)
           if (r == q) then
              ! case 1: r = q
              ! simultaneously bidiagonalize x11 and x21
              call la_qorbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_qlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_qorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_qlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_qorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 v1t(1,1) = one
                 do j = 2,q
                    v1t(1,j) = zero
                    v1t(j,1) = zero
                 end do
                 call la_qlacpy('U',q - 1,q - 1,x21(1,2),ldx21,v1t(2,2),ldv1t)
                 call la_qorglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglq,childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_qbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,work(iphi),u1, &
              ldu1,u2,ldu2,v1t,ldv1t,dum2,1,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place zero submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_qlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == p) then
              ! case 2: r = p
              ! simultaneously bidiagonalize x11 and x21
              call la_qorbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 u1(1,1) = one
                 do j = 2,p
                    u1(1,j) = zero
                    u1(j,1) = zero
                 end do
                 call la_qlacpy('L',p - 1,p - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_qorgqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,work(itaup1),work(iorgqr), &
                           lorgqr,childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_qlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_qorgqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_qlacpy('U',p,q,x11,ldx11,v1t,ldv1t)
                 call la_qorglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_qbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,work(iphi),v1t, &
              ldv1t,dum2,1,u1,ldu1,u2,ldu2,work(ib11d),work(ib11e),work(ib12d),work(ib12e) &
              ,work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd,childinfo &
                        )
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_qlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == m - p) then
              ! case 3: r = m-p
              ! simultaneously bidiagonalize x11 and x21
              call la_qorbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
                        ,work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_qlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_qorgqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 u2(1,1) = one
                 do j = 2,m - p
                    u2(1,j) = zero
                    u2(j,1) = zero
                 end do
                 call la_qlacpy('L',m - p - 1,m - p - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_qorgqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,work(itaup2),work(iorgqr) &
                           ,lorgqr,childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_qlacpy('U',m - p,q,x21,ldx21,v1t,ldv1t)
                 call la_qorglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_qbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,work(iphi), &
              dum2,1,v1t,ldv1t,u2,ldu2,u1,ldu1,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > r) then
                 do i = 1,r
                    iwork(i) = q - r + i
                 end do
                 do i = r + 1,q
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_qlapmt(.false.,p,q,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_qlapmr(.false.,q,q,v1t,ldv1t,iwork)
                 end if
              end if
           else
              ! case 4: r = m-q
              ! simultaneously bidiagonalize x11 and x21
              call la_qorbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,work(iphi),work(itaup1) &
              ,work(itaup2),work(itauq1),work(iorbdb),work(iorbdb + m),lorbdb - m,childinfo)

              ! accumulate householder reflectors
              if (wantu2 .and. m - p > 0) then
                 call la_qcopy(m - p,work(iorbdb + p),1,u2,1)
              end if
              if (wantu1 .and. p > 0) then
                 call la_qcopy(p,work(iorbdb),1,u1,1)
                 do j = 2,p
                    u1(1,j) = zero
                 end do
                 call la_qlacpy('L',p - 1,m - q - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_qorgqr(p,p,m - q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 do j = 2,m - p
                    u2(1,j) = zero
                 end do
                 call la_qlacpy('L',m - p - 1,m - q - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_qorgqr(m - p,m - p,m - q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_qlacpy('U',m - q,q,x21,ldx21,v1t,ldv1t)
                 call la_qlacpy('U',p - (m - q),q - (m - q),x11(m - q + 1,m - q + 1),ldx11,v1t(m - q + 1,m - q + &
                           1),ldv1t)
                 call la_qlacpy('U',-p + q,q - p,x21(m - q + 1,p + 1),ldx21,v1t(p + 1,p + 1),ldv1t)

                 call la_qorglq(q,q,q,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_qbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,work(iphi), &
              u2,ldu2,u1,ldu1,dum2,1,v1t,ldv1t,work(ib11d),work(ib11e),work(ib12d),work( &
              ib12e),work(ib21d),work(ib21e),work(ib22d),work(ib22e),work(ibbcsd),lbbcsd, &
                        childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (p > r) then
                 do i = 1,r
                    iwork(i) = p - r + i
                 end do
                 do i = r + 1,p
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_qlapmt(.false.,p,p,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_qlapmr(.false.,p,q,v1t,ldv1t,iwork)
                 end if
              end if
           end if
           return
     end subroutine la_qorcsd2by1

     !> CLAPMR: rearranges the rows of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.

     pure subroutine la_clapmr(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           complex(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,in,j,jj
           complex(sp) :: temp
           ! Executable Statements
           if (m <= 1) return
           do i = 1,m
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,m
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do jj = 1,n
                    temp = x(j,jj)
                    x(j,jj) = x(in,jj)
                    x(in,jj) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,m
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do jj = 1,n
                    temp = x(i,jj)
                    x(i,jj) = x(j,jj)
                    x(j,jj) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_clapmr
     !> ZLAPMR: rearranges the rows of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.

     pure subroutine la_zlapmr(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           complex(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,in,j,jj
           complex(dp) :: temp
           ! Executable Statements
           if (m <= 1) return
           do i = 1,m
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,m
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do jj = 1,n
                    temp = x(j,jj)
                    x(j,jj) = x(in,jj)
                    x(in,jj) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,m
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do jj = 1,n
                    temp = x(i,jj)
                    x(i,jj) = x(j,jj)
                    x(j,jj) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_zlapmr
     !> WLAPMR: rearranges the rows of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.

     pure subroutine la_wlapmr(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           complex(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,in,j,jj
           complex(qp) :: temp
           ! Executable Statements
           if (m <= 1) return
           do i = 1,m
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,m
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do jj = 1,n
                    temp = x(j,jj)
                    x(j,jj) = x(in,jj)
                    x(in,jj) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,m
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do jj = 1,n
                    temp = x(i,jj)
                    x(i,jj) = x(j,jj)
                    x(j,jj) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_wlapmr

     !> CLAPMT: rearranges the columns of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.

     pure subroutine la_clapmt(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           complex(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ii,j,in
           complex(sp) :: temp
           ! Executable Statements
           if (n <= 1) return
           do i = 1,n
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,n
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do ii = 1,m
                    temp = x(ii,j)
                    x(ii,j) = x(ii,in)
                    x(ii,in) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,n
                 if (k(i) > 0) go to 100
                 k(i) = -k(i)
                 j = k(i)
                 80 continue
                 if (j == i) go to 100
                 do ii = 1,m
                    temp = x(ii,i)
                    x(ii,i) = x(ii,j)
                    x(ii,j) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 80
                 100 continue
              end do
           end if
           return
     end subroutine la_clapmt
     !> ZLAPMT: rearranges the columns of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.

     pure subroutine la_zlapmt(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           complex(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ii,in,j
           complex(dp) :: temp
           ! Executable Statements
           if (n <= 1) return
           do i = 1,n
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,n
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do ii = 1,m
                    temp = x(ii,j)
                    x(ii,j) = x(ii,in)
                    x(ii,in) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,n
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do ii = 1,m
                    temp = x(ii,i)
                    x(ii,i) = x(ii,j)
                    x(ii,j) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_zlapmt
     !> WLAPMT: rearranges the columns of the M by N matrix X as specified
     !> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
     !> If FORWRD = .TRUE.,  forward permutation:
     !> X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
     !> If FORWRD = .FALSE., backward permutation:
     !> X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.

     pure subroutine la_wlapmt(forwrd,m,n,x,ldx,k)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: forwrd
           integer(ilp),intent(in) :: ldx,m,n
           ! Array Arguments
           integer(ilp),intent(inout) :: k(*)
           complex(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ii,in,j
           complex(qp) :: temp
           ! Executable Statements
           if (n <= 1) return
           do i = 1,n
              k(i) = -k(i)
           end do
           if (forwrd) then
              ! forward permutation
              do i = 1,n
                 if (k(i) > 0) go to 40
                 j = i
                 k(j) = -k(j)
                 in = k(j)
                 20 continue
                 if (k(in) > 0) go to 40
                 do ii = 1,m
                    temp = x(ii,j)
                    x(ii,j) = x(ii,in)
                    x(ii,in) = temp
                 end do
                 k(in) = -k(in)
                 j = in
                 in = k(in)
                 go to 20
                 40 continue
              end do
           else
              ! backward permutation
              do i = 1,n
                 if (k(i) > 0) go to 80
                 k(i) = -k(i)
                 j = k(i)
                 60 continue
                 if (j == i) go to 80
                 do ii = 1,m
                    temp = x(ii,i)
                    x(ii,i) = x(ii,j)
                    x(ii,j) = temp
                 end do
                 k(j) = -k(j)
                 j = k(j)
                 go to 60
                 80 continue
              end do
           end if
           return
     end subroutine la_wlapmt

     !> CUNBDB: simultaneously bidiagonalizes the blocks of an M-by-M
     !> partitioned unitary matrix X:
     !> [ B11 | B12 0  0 ]
     !> [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**H
     !> X = [-----------] = [---------] [----------------] [---------]   .
     !> [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
     !> [  0  |  0  0  I ]
     !> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     !> not the case, then X must be transposed and/or permuted. This can be
     !> done in constant time using the TRANS and SIGNS options. See CUNCSD
     !> for details.)
     !> The unitary matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     !> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     !> represented implicitly by Householder vectors.
     !> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_cunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
               ldx22,theta,phi,taup1,taup2,tauq1,tauq2,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldx11,ldx12,ldx21,ldx22,lwork,m,p,q
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           complex(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),tauq2(*),work(*)
           complex(sp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ====================================================================
           ! Parameters
           real(sp),parameter :: realone = 1.0_sp

           ! Local Scalars
           logical(lk) :: colmajor,lquery
           integer(ilp) :: i,lworkmin,lworkopt
           real(sp) :: z1,z2,z3,z4
           ! Intrinsic Functions
           intrinsic :: atan2,cos,max,min,sin
           intrinsic :: cmplx,conjg
           ! Executable Statements
           ! test input arguments
           info = 0
           colmajor = .not. la_lsame(trans,'T')
           if (.not. la_lsame(signs,'O')) then
              z1 = realone
              z2 = realone
              z3 = realone
              z4 = realone
           else
              z1 = realone
              z2 = -realone
              z3 = realone
              z4 = -realone
           end if
           lquery = lwork == -1
           if (m < 0) then
              info = -3
           else if (p < 0 .or. p > m) then
              info = -4
           else if (q < 0 .or. q > p .or. q > m - p .or. q > m - q) then
              info = -5
           else if (colmajor .and. ldx11 < max(1,p)) then
              info = -7
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
              info = -7
           else if (colmajor .and. ldx12 < max(1,p)) then
              info = -9
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
              info = -9
           else if (colmajor .and. ldx21 < max(1,m - p)) then
              info = -11
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
              info = -11
           else if (colmajor .and. ldx22 < max(1,m - p)) then
              info = -13
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
              info = -13
           end if
           ! compute workspace
           if (info == 0) then
              lworkopt = m - q
              lworkmin = m - q
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORBDB',-info)
              return
           else if (lquery) then
              return
           end if
           ! handle column-major and row-major separately
           if (colmajor) then
              ! reduce columns 1, ..., q of x11, x12, x21, and x22
              do i = 1,q
                 if (i == 1) then
                    call la_cscal(p - i + 1,cmplx(z1,0.0_sp,KIND=sp),x11(i,i),1)
                 else
                    call la_cscal(p - i + 1,cmplx(z1*cos(phi(i - 1)),0.0_sp,KIND=sp),x11(i,i), &
                              1)
                    call la_caxpy(p - i + 1,cmplx(-z1*z3*z4*sin(phi(i - 1)),0.0_sp,KIND=sp),x12( &
                              i,i - 1),1,x11(i,i),1)
                 end if
                 if (i == 1) then
                    call la_cscal(m - p - i + 1,cmplx(z2,0.0_sp,KIND=sp),x21(i,i),1)
                 else
                    call la_cscal(m - p - i + 1,cmplx(z2*cos(phi(i - 1)),0.0_sp,KIND=sp),x21(i,i), &
                               1)
                    call la_caxpy(m - p - i + 1,cmplx(-z2*z3*z4*sin(phi(i - 1)),0.0_sp,KIND=sp), &
                              x22(i,i - 1),1,x21(i,i),1)
                 end if
                 theta(i) = atan2(la_scnrm2(m - p - i + 1,x21(i,i),1),la_scnrm2(p - i + 1, &
                           x11(i,i),1))
                 if (p > i) then
                    call la_clarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
                 else if (p == i) then
                    call la_clarfgp(p - i + 1,x11(i,i),x11(i,i),1,taup1(i))
                 end if
                 x11(i,i) = cone
                 if (m - p > i) then
                    call la_clarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
                 else if (m - p == i) then
                    call la_clarfgp(m - p - i + 1,x21(i,i),x21(i,i),1,taup2(i))
                 end if
                 x21(i,i) = cone
                 if (q > i) then
                    call la_clarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1), &
                              ldx11,work)
                    call la_clarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                               ldx21,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_clarf('L',p - i + 1,m - q - i + 1,x11(i,i),1,conjg(taup1(i)),x12(i,i), &
                               ldx12,work)
                    call la_clarf('L',m - p - i + 1,m - q - i + 1,x21(i,i),1,conjg(taup2(i)),x22(i, &
                              i),ldx22,work)
                 end if
                 if (i < q) then
                    call la_cscal(q - i,cmplx(-z1*z3*sin(theta(i)),0.0_sp,KIND=sp),x11(i,i + &
                              1),ldx11)
                    call la_caxpy(q - i,cmplx(z2*z3*cos(theta(i)),0.0_sp,KIND=sp),x21(i,i + 1) &
                              ,ldx21,x11(i,i + 1),ldx11)
                 end if
                 call la_cscal(m - q - i + 1,cmplx(-z1*z4*sin(theta(i)),0.0_sp,KIND=sp),x12(i,i) &
                           ,ldx12)
                 call la_caxpy(m - q - i + 1,cmplx(z2*z4*cos(theta(i)),0.0_sp,KIND=sp),x22(i,i), &
                            ldx22,x12(i,i),ldx12)
                 if (i < q) phi(i) = atan2(la_scnrm2(q - i,x11(i,i + 1),ldx11),la_scnrm2( &
                            m - q - i + 1,x12(i,i),ldx12))
                 if (i < q) then
                    call la_clacgv(q - i,x11(i,i + 1),ldx11)
                    if (i == q - 1) then
                       call la_clarfgp(q - i,x11(i,i + 1),x11(i,i + 1),ldx11,tauq1(i))
                    else
                       call la_clarfgp(q - i,x11(i,i + 1),x11(i,i + 2),ldx11,tauq1(i))
                    end if
                    x11(i,i + 1) = cone
                 end if
                 if (m - q + 1 > i) then
                    call la_clacgv(m - q - i + 1,x12(i,i),ldx12)
                    if (m - q == i) then
                       call la_clarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                    else
                       call la_clarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))

                    end if
                 end if
                 x12(i,i) = cone
                 if (i < q) then
                    call la_clarf('R',p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x11(i + 1,i + 1), &
                              ldx11,work)
                    call la_clarf('R',m - p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 if (p > i) then
                    call la_clarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p > i) then
                    call la_clarf('R',m - p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x22(i + 1,i), &
                              ldx22,work)
                 end if
                 if (i < q) call la_clacgv(q - i,x11(i,i + 1),ldx11)
                 call la_clacgv(m - q - i + 1,x12(i,i),ldx12)
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_cscal(m - q - i + 1,cmplx(-z1*z4,0.0_sp,KIND=sp),x12(i,i),ldx12)

                 call la_clacgv(m - q - i + 1,x12(i,i),ldx12)
                 if (i >= m - q) then
                    call la_clarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                 else
                    call la_clarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))
                 end if
                 x12(i,i) = cone
                 if (p > i) then
                    call la_clarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_clarf('R',m - p - q,m - q - i + 1,x12(i,i),ldx12,tauq2(i), &
                            x22(q + 1,i),ldx22,work)
                 call la_clacgv(m - q - i + 1,x12(i,i),ldx12)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_cscal(m - p - q - i + 1,cmplx(z2*z4,0.0_sp,KIND=sp),x22(q + i,p + i),ldx22)

                 call la_clacgv(m - p - q - i + 1,x22(q + i,p + i),ldx22)
                 call la_clarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i + 1),ldx22,tauq2(p + i))

                 x22(q + i,p + i) = cone
                 call la_clarf('R',m - p - q - i,m - p - q - i + 1,x22(q + i,p + i),ldx22,tauq2(p + i),x22( &
                           q + i + 1,p + i),ldx22,work)
                 call la_clacgv(m - p - q - i + 1,x22(q + i,p + i),ldx22)
              end do
           else
              ! reduce columns 1, ..., q of x11, x12, x21, x22
              do i = 1,q
                 if (i == 1) then
                    call la_cscal(p - i + 1,cmplx(z1,0.0_sp,KIND=sp),x11(i,i),ldx11)
                 else
                    call la_cscal(p - i + 1,cmplx(z1*cos(phi(i - 1)),0.0_sp,KIND=sp),x11(i,i), &
                              ldx11)
                    call la_caxpy(p - i + 1,cmplx(-z1*z3*z4*sin(phi(i - 1)),0.0_sp,KIND=sp),x12( &
                              i - 1,i),ldx12,x11(i,i),ldx11)
                 end if
                 if (i == 1) then
                    call la_cscal(m - p - i + 1,cmplx(z2,0.0_sp,KIND=sp),x21(i,i),ldx21)

                 else
                    call la_cscal(m - p - i + 1,cmplx(z2*cos(phi(i - 1)),0.0_sp,KIND=sp),x21(i,i), &
                               ldx21)
                    call la_caxpy(m - p - i + 1,cmplx(-z2*z3*z4*sin(phi(i - 1)),0.0_sp,KIND=sp), &
                              x22(i - 1,i),ldx22,x21(i,i),ldx21)
                 end if
                 theta(i) = atan2(la_scnrm2(m - p - i + 1,x21(i,i),ldx21),la_scnrm2(p - i + 1, &
                            x11(i,i),ldx11))
                 call la_clacgv(p - i + 1,x11(i,i),ldx11)
                 call la_clacgv(m - p - i + 1,x21(i,i),ldx21)
                 call la_clarfgp(p - i + 1,x11(i,i),x11(i,i + 1),ldx11,taup1(i))
                 x11(i,i) = cone
                 if (i == m - p) then
                    call la_clarfgp(m - p - i + 1,x21(i,i),x21(i,i),ldx21,taup2(i))
                 else
                    call la_clarfgp(m - p - i + 1,x21(i,i),x21(i,i + 1),ldx21,taup2(i))
                 end if
                 x21(i,i) = cone
                 call la_clarf('R',q - i,p - i + 1,x11(i,i),ldx11,taup1(i),x11(i + 1,i),ldx11, &
                           work)
                 call la_clarf('R',m - q - i + 1,p - i + 1,x11(i,i),ldx11,taup1(i),x12(i,i), &
                           ldx12,work)
                 call la_clarf('R',q - i,m - p - i + 1,x21(i,i),ldx21,taup2(i),x21(i + 1,i), &
                           ldx21,work)
                 call la_clarf('R',m - q - i + 1,m - p - i + 1,x21(i,i),ldx21,taup2(i),x22(i,i), &
                           ldx22,work)
                 call la_clacgv(p - i + 1,x11(i,i),ldx11)
                 call la_clacgv(m - p - i + 1,x21(i,i),ldx21)
                 if (i < q) then
                    call la_cscal(q - i,cmplx(-z1*z3*sin(theta(i)),0.0_sp,KIND=sp),x11(i + 1, &
                              i),1)
                    call la_caxpy(q - i,cmplx(z2*z3*cos(theta(i)),0.0_sp,KIND=sp),x21(i + 1,i) &
                              ,1,x11(i + 1,i),1)
                 end if
                 call la_cscal(m - q - i + 1,cmplx(-z1*z4*sin(theta(i)),0.0_sp,KIND=sp),x12(i,i) &
                           ,1)
                 call la_caxpy(m - q - i + 1,cmplx(z2*z4*cos(theta(i)),0.0_sp,KIND=sp),x22(i,i), &
                            1,x12(i,i),1)
                 if (i < q) phi(i) = atan2(la_scnrm2(q - i,x11(i + 1,i),1),la_scnrm2(m - &
                           q - i + 1,x12(i,i),1))
                 if (i < q) then
                    call la_clarfgp(q - i,x11(i + 1,i),x11(i + 2,i),1,tauq1(i))
                    x11(i + 1,i) = cone
                 end if
                 call la_clarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = cone
                 if (i < q) then
                    call la_clarf('L',q - i,p - i,x11(i + 1,i),1,conjg(tauq1(i)),x11(i + 1,i + 1), &
                               ldx11,work)
                    call la_clarf('L',q - i,m - p - i,x11(i + 1,i),1,conjg(tauq1(i)),x21(i + 1,i + &
                              1),ldx21,work)
                 end if
                 call la_clarf('L',m - q - i + 1,p - i,x12(i,i),1,conjg(tauq2(i)),x12(i,i + 1), &
                           ldx12,work)
                 if (m - p > i) then
                    call la_clarf('L',m - q - i + 1,m - p - i,x12(i,i),1,conjg(tauq2(i)),x22(i,i + &
                              1),ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_cscal(m - q - i + 1,cmplx(-z1*z4,0.0_sp,KIND=sp),x12(i,i),1)
                 call la_clarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = cone
                 if (p > i) then
                    call la_clarf('L',m - q - i + 1,p - i,x12(i,i),1,conjg(tauq2(i)),x12(i,i + 1), &
                               ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_clarf('L',m - q - i + 1,m - p - q,x12(i,i),1,conjg(tauq2( &
                           i)),x22(i,q + 1),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_cscal(m - p - q - i + 1,cmplx(z2*z4,0.0_sp,KIND=sp),x22(p + i,q + i),1)

                 call la_clarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i + 1,q + i),1,tauq2(p + i))

                 x22(p + i,q + i) = cone
                 if (m - p - q /= i) then
                    call la_clarf('L',m - p - q - i + 1,m - p - q - i,x22(p + i,q + i),1,conjg(tauq2(p + i)), &
                               x22(p + i,q + i + 1),ldx22,work)
                 end if
              end do
           end if
           return
     end subroutine la_cunbdb
     !> ZUNBDB: simultaneously bidiagonalizes the blocks of an M-by-M
     !> partitioned unitary matrix X:
     !> [ B11 | B12 0  0 ]
     !> [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**H
     !> X = [-----------] = [---------] [----------------] [---------]   .
     !> [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
     !> [  0  |  0  0  I ]
     !> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     !> not the case, then X must be transposed and/or permuted. This can be
     !> done in constant time using the TRANS and SIGNS options. See ZUNCSD
     !> for details.)
     !> The unitary matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     !> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     !> represented implicitly by Householder vectors.
     !> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_zunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
               ldx22,theta,phi,taup1,taup2,tauq1,tauq2,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldx11,ldx12,ldx21,ldx22,lwork,m,p,q
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           complex(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),tauq2(*),work(*)
           complex(dp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ====================================================================
           ! Parameters
           real(dp),parameter :: realone = 1.0_dp

           ! Local Scalars
           logical(lk) :: colmajor,lquery
           integer(ilp) :: i,lworkmin,lworkopt
           real(dp) :: z1,z2,z3,z4
           ! Intrinsic Functions
           intrinsic :: atan2,cos,max,min,sin
           intrinsic :: cmplx,conjg
           ! Executable Statements
           ! test input arguments
           info = 0
           colmajor = .not. la_lsame(trans,'T')
           if (.not. la_lsame(signs,'O')) then
              z1 = realone
              z2 = realone
              z3 = realone
              z4 = realone
           else
              z1 = realone
              z2 = -realone
              z3 = realone
              z4 = -realone
           end if
           lquery = lwork == -1
           if (m < 0) then
              info = -3
           else if (p < 0 .or. p > m) then
              info = -4
           else if (q < 0 .or. q > p .or. q > m - p .or. q > m - q) then
              info = -5
           else if (colmajor .and. ldx11 < max(1,p)) then
              info = -7
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
              info = -7
           else if (colmajor .and. ldx12 < max(1,p)) then
              info = -9
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
              info = -9
           else if (colmajor .and. ldx21 < max(1,m - p)) then
              info = -11
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
              info = -11
           else if (colmajor .and. ldx22 < max(1,m - p)) then
              info = -13
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
              info = -13
           end if
           ! compute workspace
           if (info == 0) then
              lworkopt = m - q
              lworkmin = m - q
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORBDB',-info)
              return
           else if (lquery) then
              return
           end if
           ! handle column-major and row-major separately
           if (colmajor) then
              ! reduce columns 1, ..., q of x11, x12, x21, and x22
              do i = 1,q
                 if (i == 1) then
                    call la_zscal(p - i + 1,cmplx(z1,0.0_dp,KIND=dp),x11(i,i),1)
                 else
                    call la_zscal(p - i + 1,cmplx(z1*cos(phi(i - 1)),0.0_dp,KIND=dp),x11(i,i), &
                              1)
                    call la_zaxpy(p - i + 1,cmplx(-z1*z3*z4*sin(phi(i - 1)),0.0_dp,KIND=dp),x12( &
                              i,i - 1),1,x11(i,i),1)
                 end if
                 if (i == 1) then
                    call la_zscal(m - p - i + 1,cmplx(z2,0.0_dp,KIND=dp),x21(i,i),1)
                 else
                    call la_zscal(m - p - i + 1,cmplx(z2*cos(phi(i - 1)),0.0_dp,KIND=dp),x21(i,i), &
                               1)
                    call la_zaxpy(m - p - i + 1,cmplx(-z2*z3*z4*sin(phi(i - 1)),0.0_dp,KIND=dp), &
                              x22(i,i - 1),1,x21(i,i),1)
                 end if
                 theta(i) = atan2(la_dznrm2(m - p - i + 1,x21(i,i),1),la_dznrm2(p - i + 1, &
                           x11(i,i),1))
                 if (p > i) then
                    call la_zlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
                 else if (p == i) then
                    call la_zlarfgp(p - i + 1,x11(i,i),x11(i,i),1,taup1(i))
                 end if
                 x11(i,i) = cone
                 if (m - p > i) then
                    call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
                 else if (m - p == i) then
                    call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i,i),1,taup2(i))
                 end if
                 x21(i,i) = cone
                 if (q > i) then
                    call la_zlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1), &
                              ldx11,work)
                    call la_zlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                               ldx21,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_zlarf('L',p - i + 1,m - q - i + 1,x11(i,i),1,conjg(taup1(i)),x12(i,i), &
                               ldx12,work)
                    call la_zlarf('L',m - p - i + 1,m - q - i + 1,x21(i,i),1,conjg(taup2(i)),x22(i, &
                              i),ldx22,work)
                 end if
                 if (i < q) then
                    call la_zscal(q - i,cmplx(-z1*z3*sin(theta(i)),0.0_dp,KIND=dp),x11(i,i + &
                              1),ldx11)
                    call la_zaxpy(q - i,cmplx(z2*z3*cos(theta(i)),0.0_dp,KIND=dp),x21(i,i + 1) &
                              ,ldx21,x11(i,i + 1),ldx11)
                 end if
                 call la_zscal(m - q - i + 1,cmplx(-z1*z4*sin(theta(i)),0.0_dp,KIND=dp),x12(i,i) &
                           ,ldx12)
                 call la_zaxpy(m - q - i + 1,cmplx(z2*z4*cos(theta(i)),0.0_dp,KIND=dp),x22(i,i), &
                            ldx22,x12(i,i),ldx12)
                 if (i < q) phi(i) = atan2(la_dznrm2(q - i,x11(i,i + 1),ldx11),la_dznrm2( &
                            m - q - i + 1,x12(i,i),ldx12))
                 if (i < q) then
                    call la_zlacgv(q - i,x11(i,i + 1),ldx11)
                    if (i == q - 1) then
                       call la_zlarfgp(q - i,x11(i,i + 1),x11(i,i + 1),ldx11,tauq1(i))
                    else
                       call la_zlarfgp(q - i,x11(i,i + 1),x11(i,i + 2),ldx11,tauq1(i))
                    end if
                    x11(i,i + 1) = cone
                 end if
                 if (m - q + 1 > i) then
                    call la_zlacgv(m - q - i + 1,x12(i,i),ldx12)
                    if (m - q == i) then
                       call la_zlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                    else
                       call la_zlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))

                    end if
                 end if
                 x12(i,i) = cone
                 if (i < q) then
                    call la_zlarf('R',p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x11(i + 1,i + 1), &
                              ldx11,work)
                    call la_zlarf('R',m - p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 if (p > i) then
                    call la_zlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p > i) then
                    call la_zlarf('R',m - p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x22(i + 1,i), &
                              ldx22,work)
                 end if
                 if (i < q) call la_zlacgv(q - i,x11(i,i + 1),ldx11)
                 call la_zlacgv(m - q - i + 1,x12(i,i),ldx12)
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_zscal(m - q - i + 1,cmplx(-z1*z4,0.0_dp,KIND=dp),x12(i,i),ldx12)

                 call la_zlacgv(m - q - i + 1,x12(i,i),ldx12)
                 if (i >= m - q) then
                    call la_zlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                 else
                    call la_zlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))
                 end if
                 x12(i,i) = cone
                 if (p > i) then
                    call la_zlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_zlarf('R',m - p - q,m - q - i + 1,x12(i,i),ldx12,tauq2(i), &
                            x22(q + 1,i),ldx22,work)
                 call la_zlacgv(m - q - i + 1,x12(i,i),ldx12)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_zscal(m - p - q - i + 1,cmplx(z2*z4,0.0_dp,KIND=dp),x22(q + i,p + i),ldx22)

                 call la_zlacgv(m - p - q - i + 1,x22(q + i,p + i),ldx22)
                 call la_zlarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i + 1),ldx22,tauq2(p + i))

                 x22(q + i,p + i) = cone
                 call la_zlarf('R',m - p - q - i,m - p - q - i + 1,x22(q + i,p + i),ldx22,tauq2(p + i),x22( &
                           q + i + 1,p + i),ldx22,work)
                 call la_zlacgv(m - p - q - i + 1,x22(q + i,p + i),ldx22)
              end do
           else
              ! reduce columns 1, ..., q of x11, x12, x21, x22
              do i = 1,q
                 if (i == 1) then
                    call la_zscal(p - i + 1,cmplx(z1,0.0_dp,KIND=dp),x11(i,i),ldx11)
                 else
                    call la_zscal(p - i + 1,cmplx(z1*cos(phi(i - 1)),0.0_dp,KIND=dp),x11(i,i), &
                              ldx11)
                    call la_zaxpy(p - i + 1,cmplx(-z1*z3*z4*sin(phi(i - 1)),0.0_dp,KIND=dp),x12( &
                              i - 1,i),ldx12,x11(i,i),ldx11)
                 end if
                 if (i == 1) then
                    call la_zscal(m - p - i + 1,cmplx(z2,0.0_dp,KIND=dp),x21(i,i),ldx21)

                 else
                    call la_zscal(m - p - i + 1,cmplx(z2*cos(phi(i - 1)),0.0_dp,KIND=dp),x21(i,i), &
                               ldx21)
                    call la_zaxpy(m - p - i + 1,cmplx(-z2*z3*z4*sin(phi(i - 1)),0.0_dp,KIND=dp), &
                              x22(i - 1,i),ldx22,x21(i,i),ldx21)
                 end if
                 theta(i) = atan2(la_dznrm2(m - p - i + 1,x21(i,i),ldx21),la_dznrm2(p - i + 1, &
                            x11(i,i),ldx11))
                 call la_zlacgv(p - i + 1,x11(i,i),ldx11)
                 call la_zlacgv(m - p - i + 1,x21(i,i),ldx21)
                 call la_zlarfgp(p - i + 1,x11(i,i),x11(i,i + 1),ldx11,taup1(i))
                 x11(i,i) = cone
                 if (i == m - p) then
                    call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i,i),ldx21,taup2(i))
                 else
                    call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i,i + 1),ldx21,taup2(i))
                 end if
                 x21(i,i) = cone
                 call la_zlarf('R',q - i,p - i + 1,x11(i,i),ldx11,taup1(i),x11(i + 1,i),ldx11, &
                           work)
                 call la_zlarf('R',m - q - i + 1,p - i + 1,x11(i,i),ldx11,taup1(i),x12(i,i), &
                           ldx12,work)
                 call la_zlarf('R',q - i,m - p - i + 1,x21(i,i),ldx21,taup2(i),x21(i + 1,i), &
                           ldx21,work)
                 call la_zlarf('R',m - q - i + 1,m - p - i + 1,x21(i,i),ldx21,taup2(i),x22(i,i), &
                           ldx22,work)
                 call la_zlacgv(p - i + 1,x11(i,i),ldx11)
                 call la_zlacgv(m - p - i + 1,x21(i,i),ldx21)
                 if (i < q) then
                    call la_zscal(q - i,cmplx(-z1*z3*sin(theta(i)),0.0_dp,KIND=dp),x11(i + 1, &
                              i),1)
                    call la_zaxpy(q - i,cmplx(z2*z3*cos(theta(i)),0.0_dp,KIND=dp),x21(i + 1,i) &
                              ,1,x11(i + 1,i),1)
                 end if
                 call la_zscal(m - q - i + 1,cmplx(-z1*z4*sin(theta(i)),0.0_dp,KIND=dp),x12(i,i) &
                           ,1)
                 call la_zaxpy(m - q - i + 1,cmplx(z2*z4*cos(theta(i)),0.0_dp,KIND=dp),x22(i,i), &
                            1,x12(i,i),1)
                 if (i < q) phi(i) = atan2(la_dznrm2(q - i,x11(i + 1,i),1),la_dznrm2(m - &
                           q - i + 1,x12(i,i),1))
                 if (i < q) then
                    call la_zlarfgp(q - i,x11(i + 1,i),x11(i + 2,i),1,tauq1(i))
                    x11(i + 1,i) = cone
                 end if
                 call la_zlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = cone
                 if (i < q) then
                    call la_zlarf('L',q - i,p - i,x11(i + 1,i),1,conjg(tauq1(i)),x11(i + 1,i + 1), &
                               ldx11,work)
                    call la_zlarf('L',q - i,m - p - i,x11(i + 1,i),1,conjg(tauq1(i)),x21(i + 1,i + &
                              1),ldx21,work)
                 end if
                 call la_zlarf('L',m - q - i + 1,p - i,x12(i,i),1,conjg(tauq2(i)),x12(i,i + 1), &
                           ldx12,work)
                 if (m - p > i) then
                    call la_zlarf('L',m - q - i + 1,m - p - i,x12(i,i),1,conjg(tauq2(i)),x22(i,i + &
                              1),ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_zscal(m - q - i + 1,cmplx(-z1*z4,0.0_dp,KIND=dp),x12(i,i),1)
                 call la_zlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = cone
                 if (p > i) then
                    call la_zlarf('L',m - q - i + 1,p - i,x12(i,i),1,conjg(tauq2(i)),x12(i,i + 1), &
                               ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_zlarf('L',m - q - i + 1,m - p - q,x12(i,i),1,conjg(tauq2( &
                           i)),x22(i,q + 1),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_zscal(m - p - q - i + 1,cmplx(z2*z4,0.0_dp,KIND=dp),x22(p + i,q + i),1)

                 call la_zlarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i + 1,q + i),1,tauq2(p + i))

                 x22(p + i,q + i) = cone
                 if (m - p - q /= i) then
                    call la_zlarf('L',m - p - q - i + 1,m - p - q - i,x22(p + i,q + i),1,conjg(tauq2(p + i)), &
                               x22(p + i,q + i + 1),ldx22,work)
                 end if
              end do
           end if
           return
     end subroutine la_zunbdb
     !> WUNBDB: simultaneously bidiagonalizes the blocks of an M-by-M
     !> partitioned unitary matrix X:
     !> [ B11 | B12 0  0 ]
     !> [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**H
     !> X = [-----------] = [---------] [----------------] [---------]   .
     !> [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
     !> [  0  |  0  0  I ]
     !> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
     !> not the case, then X must be transposed and/or permuted. This can be
     !> done in constant time using the TRANS and SIGNS options. See WUNCSD
     !> for details.)
     !> The unitary matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
     !> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
     !> represented implicitly by Householder vectors.
     !> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_wunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
               ldx22,theta,phi,taup1,taup2,tauq1,tauq2,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldx11,ldx12,ldx21,ldx22,lwork,m,p,q
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           complex(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),tauq2(*),work(*)
           complex(qp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ====================================================================
           ! Parameters
           real(qp),parameter :: realone = 1.0_qp

           ! Local Scalars
           logical(lk) :: colmajor,lquery
           integer(ilp) :: i,lworkmin,lworkopt
           real(qp) :: z1,z2,z3,z4
           ! Intrinsic Functions
           intrinsic :: atan2,cos,max,min,sin
           intrinsic :: cmplx,conjg
           ! Executable Statements
           ! test input arguments
           info = 0
           colmajor = .not. la_lsame(trans,'T')
           if (.not. la_lsame(signs,'O')) then
              z1 = realone
              z2 = realone
              z3 = realone
              z4 = realone
           else
              z1 = realone
              z2 = -realone
              z3 = realone
              z4 = -realone
           end if
           lquery = lwork == -1
           if (m < 0) then
              info = -3
           else if (p < 0 .or. p > m) then
              info = -4
           else if (q < 0 .or. q > p .or. q > m - p .or. q > m - q) then
              info = -5
           else if (colmajor .and. ldx11 < max(1,p)) then
              info = -7
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
              info = -7
           else if (colmajor .and. ldx12 < max(1,p)) then
              info = -9
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
              info = -9
           else if (colmajor .and. ldx21 < max(1,m - p)) then
              info = -11
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
              info = -11
           else if (colmajor .and. ldx22 < max(1,m - p)) then
              info = -13
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
              info = -13
           end if
           ! compute workspace
           if (info == 0) then
              lworkopt = m - q
              lworkmin = m - q
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORBDB',-info)
              return
           else if (lquery) then
              return
           end if
           ! handle column-major and row-major separately
           if (colmajor) then
              ! reduce columns 1, ..., q of x11, x12, x21, and x22
              do i = 1,q
                 if (i == 1) then
                    call la_wscal(p - i + 1,cmplx(z1,0.0_qp,KIND=qp),x11(i,i),1)
                 else
                    call la_wscal(p - i + 1,cmplx(z1*cos(phi(i - 1)),0.0_qp,KIND=qp),x11(i,i), &
                              1)
                    call la_waxpy(p - i + 1,cmplx(-z1*z3*z4*sin(phi(i - 1)),0.0_qp,KIND=qp),x12( &
                              i,i - 1),1,x11(i,i),1)
                 end if
                 if (i == 1) then
                    call la_wscal(m - p - i + 1,cmplx(z2,0.0_qp,KIND=qp),x21(i,i),1)
                 else
                    call la_wscal(m - p - i + 1,cmplx(z2*cos(phi(i - 1)),0.0_qp,KIND=qp),x21(i,i), &
                               1)
                    call la_waxpy(m - p - i + 1,cmplx(-z2*z3*z4*sin(phi(i - 1)),0.0_qp,KIND=qp), &
                              x22(i,i - 1),1,x21(i,i),1)
                 end if
                 theta(i) = atan2(la_qwnrm2(m - p - i + 1,x21(i,i),1),la_qwnrm2(p - i + 1, &
                           x11(i,i),1))
                 if (p > i) then
                    call la_wlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
                 else if (p == i) then
                    call la_wlarfgp(p - i + 1,x11(i,i),x11(i,i),1,taup1(i))
                 end if
                 x11(i,i) = cone
                 if (m - p > i) then
                    call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
                 else if (m - p == i) then
                    call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i,i),1,taup2(i))
                 end if
                 x21(i,i) = cone
                 if (q > i) then
                    call la_wlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1), &
                              ldx11,work)
                    call la_wlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                               ldx21,work)
                 end if
                 if (m - q + 1 > i) then
                    call la_wlarf('L',p - i + 1,m - q - i + 1,x11(i,i),1,conjg(taup1(i)),x12(i,i), &
                               ldx12,work)
                    call la_wlarf('L',m - p - i + 1,m - q - i + 1,x21(i,i),1,conjg(taup2(i)),x22(i, &
                              i),ldx22,work)
                 end if
                 if (i < q) then
                    call la_wscal(q - i,cmplx(-z1*z3*sin(theta(i)),0.0_qp,KIND=qp),x11(i,i + &
                              1),ldx11)
                    call la_waxpy(q - i,cmplx(z2*z3*cos(theta(i)),0.0_qp,KIND=qp),x21(i,i + 1) &
                              ,ldx21,x11(i,i + 1),ldx11)
                 end if
                 call la_wscal(m - q - i + 1,cmplx(-z1*z4*sin(theta(i)),0.0_qp,KIND=qp),x12(i,i) &
                           ,ldx12)
                 call la_waxpy(m - q - i + 1,cmplx(z2*z4*cos(theta(i)),0.0_qp,KIND=qp),x22(i,i), &
                            ldx22,x12(i,i),ldx12)
                 if (i < q) phi(i) = atan2(la_qwnrm2(q - i,x11(i,i + 1),ldx11),la_qwnrm2( &
                            m - q - i + 1,x12(i,i),ldx12))
                 if (i < q) then
                    call la_wlacgv(q - i,x11(i,i + 1),ldx11)
                    if (i == q - 1) then
                       call la_wlarfgp(q - i,x11(i,i + 1),x11(i,i + 1),ldx11,tauq1(i))
                    else
                       call la_wlarfgp(q - i,x11(i,i + 1),x11(i,i + 2),ldx11,tauq1(i))
                    end if
                    x11(i,i + 1) = cone
                 end if
                 if (m - q + 1 > i) then
                    call la_wlacgv(m - q - i + 1,x12(i,i),ldx12)
                    if (m - q == i) then
                       call la_wlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                    else
                       call la_wlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))

                    end if
                 end if
                 x12(i,i) = cone
                 if (i < q) then
                    call la_wlarf('R',p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x11(i + 1,i + 1), &
                              ldx11,work)
                    call la_wlarf('R',m - p - i,q - i,x11(i,i + 1),ldx11,tauq1(i),x21(i + 1,i + 1), &
                              ldx21,work)
                 end if
                 if (p > i) then
                    call la_wlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p > i) then
                    call la_wlarf('R',m - p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x22(i + 1,i), &
                              ldx22,work)
                 end if
                 if (i < q) call la_wlacgv(q - i,x11(i,i + 1),ldx11)
                 call la_wlacgv(m - q - i + 1,x12(i,i),ldx12)
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_wscal(m - q - i + 1,cmplx(-z1*z4,0.0_qp,KIND=qp),x12(i,i),ldx12)

                 call la_wlacgv(m - q - i + 1,x12(i,i),ldx12)
                 if (i >= m - q) then
                    call la_wlarfgp(m - q - i + 1,x12(i,i),x12(i,i),ldx12,tauq2(i))
                 else
                    call la_wlarfgp(m - q - i + 1,x12(i,i),x12(i,i + 1),ldx12,tauq2(i))
                 end if
                 x12(i,i) = cone
                 if (p > i) then
                    call la_wlarf('R',p - i,m - q - i + 1,x12(i,i),ldx12,tauq2(i),x12(i + 1,i), &
                              ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_wlarf('R',m - p - q,m - q - i + 1,x12(i,i),ldx12,tauq2(i), &
                            x22(q + 1,i),ldx22,work)
                 call la_wlacgv(m - q - i + 1,x12(i,i),ldx12)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_wscal(m - p - q - i + 1,cmplx(z2*z4,0.0_qp,KIND=qp),x22(q + i,p + i),ldx22)

                 call la_wlacgv(m - p - q - i + 1,x22(q + i,p + i),ldx22)
                 call la_wlarfgp(m - p - q - i + 1,x22(q + i,p + i),x22(q + i,p + i + 1),ldx22,tauq2(p + i))

                 x22(q + i,p + i) = cone
                 call la_wlarf('R',m - p - q - i,m - p - q - i + 1,x22(q + i,p + i),ldx22,tauq2(p + i),x22( &
                           q + i + 1,p + i),ldx22,work)
                 call la_wlacgv(m - p - q - i + 1,x22(q + i,p + i),ldx22)
              end do
           else
              ! reduce columns 1, ..., q of x11, x12, x21, x22
              do i = 1,q
                 if (i == 1) then
                    call la_wscal(p - i + 1,cmplx(z1,0.0_qp,KIND=qp),x11(i,i),ldx11)
                 else
                    call la_wscal(p - i + 1,cmplx(z1*cos(phi(i - 1)),0.0_qp,KIND=qp),x11(i,i), &
                              ldx11)
                    call la_waxpy(p - i + 1,cmplx(-z1*z3*z4*sin(phi(i - 1)),0.0_qp,KIND=qp),x12( &
                              i - 1,i),ldx12,x11(i,i),ldx11)
                 end if
                 if (i == 1) then
                    call la_wscal(m - p - i + 1,cmplx(z2,0.0_qp,KIND=qp),x21(i,i),ldx21)

                 else
                    call la_wscal(m - p - i + 1,cmplx(z2*cos(phi(i - 1)),0.0_qp,KIND=qp),x21(i,i), &
                               ldx21)
                    call la_waxpy(m - p - i + 1,cmplx(-z2*z3*z4*sin(phi(i - 1)),0.0_qp,KIND=qp), &
                              x22(i - 1,i),ldx22,x21(i,i),ldx21)
                 end if
                 theta(i) = atan2(la_qwnrm2(m - p - i + 1,x21(i,i),ldx21),la_qwnrm2(p - i + 1, &
                            x11(i,i),ldx11))
                 call la_wlacgv(p - i + 1,x11(i,i),ldx11)
                 call la_wlacgv(m - p - i + 1,x21(i,i),ldx21)
                 call la_wlarfgp(p - i + 1,x11(i,i),x11(i,i + 1),ldx11,taup1(i))
                 x11(i,i) = cone
                 if (i == m - p) then
                    call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i,i),ldx21,taup2(i))
                 else
                    call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i,i + 1),ldx21,taup2(i))
                 end if
                 x21(i,i) = cone
                 call la_wlarf('R',q - i,p - i + 1,x11(i,i),ldx11,taup1(i),x11(i + 1,i),ldx11, &
                           work)
                 call la_wlarf('R',m - q - i + 1,p - i + 1,x11(i,i),ldx11,taup1(i),x12(i,i), &
                           ldx12,work)
                 call la_wlarf('R',q - i,m - p - i + 1,x21(i,i),ldx21,taup2(i),x21(i + 1,i), &
                           ldx21,work)
                 call la_wlarf('R',m - q - i + 1,m - p - i + 1,x21(i,i),ldx21,taup2(i),x22(i,i), &
                           ldx22,work)
                 call la_wlacgv(p - i + 1,x11(i,i),ldx11)
                 call la_wlacgv(m - p - i + 1,x21(i,i),ldx21)
                 if (i < q) then
                    call la_wscal(q - i,cmplx(-z1*z3*sin(theta(i)),0.0_qp,KIND=qp),x11(i + 1, &
                              i),1)
                    call la_waxpy(q - i,cmplx(z2*z3*cos(theta(i)),0.0_qp,KIND=qp),x21(i + 1,i) &
                              ,1,x11(i + 1,i),1)
                 end if
                 call la_wscal(m - q - i + 1,cmplx(-z1*z4*sin(theta(i)),0.0_qp,KIND=qp),x12(i,i) &
                           ,1)
                 call la_waxpy(m - q - i + 1,cmplx(z2*z4*cos(theta(i)),0.0_qp,KIND=qp),x22(i,i), &
                            1,x12(i,i),1)
                 if (i < q) phi(i) = atan2(la_qwnrm2(q - i,x11(i + 1,i),1),la_qwnrm2(m - &
                           q - i + 1,x12(i,i),1))
                 if (i < q) then
                    call la_wlarfgp(q - i,x11(i + 1,i),x11(i + 2,i),1,tauq1(i))
                    x11(i + 1,i) = cone
                 end if
                 call la_wlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = cone
                 if (i < q) then
                    call la_wlarf('L',q - i,p - i,x11(i + 1,i),1,conjg(tauq1(i)),x11(i + 1,i + 1), &
                               ldx11,work)
                    call la_wlarf('L',q - i,m - p - i,x11(i + 1,i),1,conjg(tauq1(i)),x21(i + 1,i + &
                              1),ldx21,work)
                 end if
                 call la_wlarf('L',m - q - i + 1,p - i,x12(i,i),1,conjg(tauq2(i)),x12(i,i + 1), &
                           ldx12,work)
                 if (m - p > i) then
                    call la_wlarf('L',m - q - i + 1,m - p - i,x12(i,i),1,conjg(tauq2(i)),x22(i,i + &
                              1),ldx22,work)
                 end if
              end do
              ! reduce columns q + 1, ..., p of x12, x22
              do i = q + 1,p
                 call la_wscal(m - q - i + 1,cmplx(-z1*z4,0.0_qp,KIND=qp),x12(i,i),1)
                 call la_wlarfgp(m - q - i + 1,x12(i,i),x12(i + 1,i),1,tauq2(i))
                 x12(i,i) = cone
                 if (p > i) then
                    call la_wlarf('L',m - q - i + 1,p - i,x12(i,i),1,conjg(tauq2(i)),x12(i,i + 1), &
                               ldx12,work)
                 end if
                 if (m - p - q >= 1) call la_wlarf('L',m - q - i + 1,m - p - q,x12(i,i),1,conjg(tauq2( &
                           i)),x22(i,q + 1),ldx22,work)
              end do
              ! reduce columns p + 1, ..., m - q of x12, x22
              do i = 1,m - p - q
                 call la_wscal(m - p - q - i + 1,cmplx(z2*z4,0.0_qp,KIND=qp),x22(p + i,q + i),1)

                 call la_wlarfgp(m - p - q - i + 1,x22(p + i,q + i),x22(p + i + 1,q + i),1,tauq2(p + i))

                 x22(p + i,q + i) = cone
                 if (m - p - q /= i) then
                    call la_wlarf('L',m - p - q - i + 1,m - p - q - i,x22(p + i,q + i),1,conjg(tauq2(p + i)), &
                               x22(p + i,q + i + 1),ldx22,work)
                 end if
              end do
           end if
           return
     end subroutine la_wunbdb

     !> CUNBDB6: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then the zero vector is returned.

     pure subroutine la_cunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           complex(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: alphasq = 0.01_sp
           real(sp),parameter :: realone = 1.0_sp
           real(sp),parameter :: realzero = 0.0_sp

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: normsq1,normsq2,scl1,scl2,ssq1,ssq2
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('CUNBDB6',-info)
              return
           end if
           ! first, project x onto the orthogonal complement of q's column
           ! space
           scl1 = realzero
           ssq1 = realone
           call la_classq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_classq(m2,x2,incx2,scl2,ssq2)
           normsq1 = scl1**2*ssq1 + scl2**2*ssq2
           if (m1 == 0) then
              do i = 1,n
                 work(i) = czero
              end do
           else
              call la_cgemv('C',m1,n,cone,q1,ldq1,x1,incx1,czero,work,1)
           end if
           call la_cgemv('C',m2,n,cone,q2,ldq2,x2,incx2,cone,work,1)
           call la_cgemv('N',m1,n,cnegone,q1,ldq1,work,1,cone,x1,incx1)
           call la_cgemv('N',m2,n,cnegone,q2,ldq2,work,1,cone,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_classq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_classq(m2,x2,incx2,scl2,ssq2)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if projection is sufficiently large in norm, then stop.
           ! if projection is czero, then stop.
           ! otherwise, project again.
           if (normsq2 >= alphasq*normsq1) then
              return
           end if
           if (normsq2 == czero) then
              return
           end if
           normsq1 = normsq2
           do i = 1,n
              work(i) = czero
           end do
           if (m1 == 0) then
              do i = 1,n
                 work(i) = czero
              end do
           else
              call la_cgemv('C',m1,n,cone,q1,ldq1,x1,incx1,czero,work,1)
           end if
           call la_cgemv('C',m2,n,cone,q2,ldq2,x2,incx2,cone,work,1)
           call la_cgemv('N',m1,n,cnegone,q1,ldq1,work,1,cone,x1,incx1)
           call la_cgemv('N',m2,n,cnegone,q2,ldq2,work,1,cone,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_classq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_classq(m1,x1,incx1,scl1,ssq1)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if second projection is sufficiently large in norm, then do
           ! nothing more. alternatively, if it shrunk significantly, then
           ! truncate it to czero.
           if (normsq2 < alphasq*normsq1) then
              do i = 1,m1
                 x1(i) = czero
              end do
              do i = 1,m2
                 x2(i) = czero
              end do
           end if
           return
     end subroutine la_cunbdb6
     !> ZUNBDB6: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then the zero vector is returned.

     pure subroutine la_zunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           complex(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: alphasq = 0.01_dp
           real(dp),parameter :: realone = 1.0_dp
           real(dp),parameter :: realzero = 0.0_dp

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: normsq1,normsq2,scl1,scl2,ssq1,ssq2
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('ZUNBDB6',-info)
              return
           end if
           ! first, project x onto the orthogonal complement of q's column
           ! space
           scl1 = realzero
           ssq1 = realone
           call la_zlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_zlassq(m2,x2,incx2,scl2,ssq2)
           normsq1 = scl1**2*ssq1 + scl2**2*ssq2
           if (m1 == 0) then
              do i = 1,n
                 work(i) = czero
              end do
           else
              call la_zgemv('C',m1,n,cone,q1,ldq1,x1,incx1,czero,work,1)
           end if
           call la_zgemv('C',m2,n,cone,q2,ldq2,x2,incx2,cone,work,1)
           call la_zgemv('N',m1,n,cnegone,q1,ldq1,work,1,cone,x1,incx1)
           call la_zgemv('N',m2,n,cnegone,q2,ldq2,work,1,cone,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_zlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_zlassq(m2,x2,incx2,scl2,ssq2)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if projection is sufficiently large in norm, then stop.
           ! if projection is czero, then stop.
           ! otherwise, project again.
           if (normsq2 >= alphasq*normsq1) then
              return
           end if
           if (normsq2 == czero) then
              return
           end if
           normsq1 = normsq2
           do i = 1,n
              work(i) = czero
           end do
           if (m1 == 0) then
              do i = 1,n
                 work(i) = czero
              end do
           else
              call la_zgemv('C',m1,n,cone,q1,ldq1,x1,incx1,czero,work,1)
           end if
           call la_zgemv('C',m2,n,cone,q2,ldq2,x2,incx2,cone,work,1)
           call la_zgemv('N',m1,n,cnegone,q1,ldq1,work,1,cone,x1,incx1)
           call la_zgemv('N',m2,n,cnegone,q2,ldq2,work,1,cone,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_zlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_zlassq(m1,x1,incx1,scl1,ssq1)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if second projection is sufficiently large in norm, then do
           ! nothing more. alternatively, if it shrunk significantly, then
           ! truncate it to czero.
           if (normsq2 < alphasq*normsq1) then
              do i = 1,m1
                 x1(i) = czero
              end do
              do i = 1,m2
                 x2(i) = czero
              end do
           end if
           return
     end subroutine la_zunbdb6
     !> WUNBDB6: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then the zero vector is returned.

     pure subroutine la_wunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           complex(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: alphasq = 0.01_qp
           real(qp),parameter :: realone = 1.0_qp
           real(qp),parameter :: realzero = 0.0_qp

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: normsq1,normsq2,scl1,scl2,ssq1,ssq2
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('WUNBDB6',-info)
              return
           end if
           ! first, project x onto the orthogonal complement of q's column
           ! space
           scl1 = realzero
           ssq1 = realone
           call la_wlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_wlassq(m2,x2,incx2,scl2,ssq2)
           normsq1 = scl1**2*ssq1 + scl2**2*ssq2
           if (m1 == 0) then
              do i = 1,n
                 work(i) = czero
              end do
           else
              call la_wgemv('C',m1,n,cone,q1,ldq1,x1,incx1,czero,work,1)
           end if
           call la_wgemv('C',m2,n,cone,q2,ldq2,x2,incx2,cone,work,1)
           call la_wgemv('N',m1,n,cnegone,q1,ldq1,work,1,cone,x1,incx1)
           call la_wgemv('N',m2,n,cnegone,q2,ldq2,work,1,cone,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_wlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_wlassq(m2,x2,incx2,scl2,ssq2)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if projection is sufficiently large in norm, then stop.
           ! if projection is czero, then stop.
           ! otherwise, project again.
           if (normsq2 >= alphasq*normsq1) then
              return
           end if
           if (normsq2 == czero) then
              return
           end if
           normsq1 = normsq2
           do i = 1,n
              work(i) = czero
           end do
           if (m1 == 0) then
              do i = 1,n
                 work(i) = czero
              end do
           else
              call la_wgemv('C',m1,n,cone,q1,ldq1,x1,incx1,czero,work,1)
           end if
           call la_wgemv('C',m2,n,cone,q2,ldq2,x2,incx2,cone,work,1)
           call la_wgemv('N',m1,n,cnegone,q1,ldq1,work,1,cone,x1,incx1)
           call la_wgemv('N',m2,n,cnegone,q2,ldq2,work,1,cone,x2,incx2)
           scl1 = realzero
           ssq1 = realone
           call la_wlassq(m1,x1,incx1,scl1,ssq1)
           scl2 = realzero
           ssq2 = realone
           call la_wlassq(m1,x1,incx1,scl1,ssq1)
           normsq2 = scl1**2*ssq1 + scl2**2*ssq2
           ! if second projection is sufficiently large in norm, then do
           ! nothing more. alternatively, if it shrunk significantly, then
           ! truncate it to czero.
           if (normsq2 < alphasq*normsq1) then
              do i = 1,m1
                 x1(i) = czero
              end do
              do i = 1,m2
                 x2(i) = czero
              end do
           end if
           return
     end subroutine la_wunbdb6

     !> CBBCSD: computes the CS decomposition of a unitary matrix in
     !> bidiagonal-block form,
     !> [ B11 | B12 0  0 ]
     !> [  0  |  0 -I  0 ]
     !> X = [----------------]
     !> [ B21 | B22 0  0 ]
     !> [  0  |  0  0  I ]
     !> [  C | -S  0  0 ]
     !> [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**H
     !> = [---------] [---------------] [---------]   .
     !> [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
     !> [  0 |  0  0  I ]
     !> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
     !> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
     !> transposed and/or permuted. This can be done in constant time using
     !> the TRANS and SIGNS options. See CUNCSD for details.)
     !> The bidiagonal matrices B11, B12, B21, and B22 are represented
     !> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
     !> The unitary matrices U1, U2, V1T, and V2T are input/output.
     !> The input matrices are pre- or post-multiplied by the appropriate
     !> singular vector matrices.

     pure subroutine la_cbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,phi,u1, &
     ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,b11d,b11e,b12d,b12e,b21d,b21e,b22d,b22e,rwork, &
               lrwork,info)
        use la_constants_sp,only:zero,one,ten,cnegone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,lrwork,m,p,q
           ! Array Arguments
           real(sp),intent(out) :: b11d(*),b11e(*),b12d(*),b12e(*),b21d(*),b21e(*),b22d(*), &
                      b22e(*),rwork(*)
           real(sp),intent(inout) :: phi(*),theta(*)
           complex(sp),intent(inout) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*)

        ! ===================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 6
           real(sp),parameter :: hundred = 100.0_sp
           real(sp),parameter :: meighth = -0.125_sp
           real(sp),parameter :: piover2 = 1.57079632679489661923132169163975144210_sp

           ! Local Scalars
           logical(lk) :: colmajor,lquery,restart11,restart12,restart21,restart22,wantu1, &
                     wantu2,wantv1t,wantv2t
           integer(ilp) :: i,imin,imax,iter,iu1cs,iu1sn,iu2cs,iu2sn,iv1tcs,iv1tsn, &
                     iv2tcs,iv2tsn,j,lrworkmin,lrworkopt,maxit,mini
           real(sp) :: b11bulge,b12bulge,b21bulge,b22bulge,dummy,eps,mu,nu,r,sigma11, &
                     sigma21,temp,thetamax,thetamin,thresh,tol,tolmul,unfl,x1,x2,y1,y2
           ! Intrinsic Functions
           intrinsic :: abs,atan2,cos,max,min,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lrwork == -1
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           if (m < 0) then
              info = -6
           else if (p < 0 .or. p > m) then
              info = -7
           else if (q < 0 .or. q > m) then
              info = -8
           else if (q > p .or. q > m - p .or. q > m - q) then
              info = -8
           else if (wantu1 .and. ldu1 < p) then
              info = -12
           else if (wantu2 .and. ldu2 < m - p) then
              info = -14
           else if (wantv1t .and. ldv1t < q) then
              info = -16
           else if (wantv2t .and. ldv2t < m - q) then
              info = -18
           end if
           ! quick return if q = 0
           if (info == 0 .and. q == 0) then
              lrworkmin = 1
              rwork(1) = lrworkmin
              return
           end if
           ! compute workspace
           if (info == 0) then
              iu1cs = 1
              iu1sn = iu1cs + q
              iu2cs = iu1sn + q
              iu2sn = iu2cs + q
              iv1tcs = iu2sn + q
              iv1tsn = iv1tcs + q
              iv2tcs = iv1tsn + q
              iv2tsn = iv2tcs + q
              lrworkopt = iv2tsn + q - 1
              lrworkmin = lrworkopt
              rwork(1) = lrworkopt
              if (lrwork < lrworkmin .and. .not. lquery) then
                 info = -28
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CBBCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! get machine constants
           eps = la_slamch('EPSILON')
           unfl = la_slamch('SAFE MINIMUM')
           tolmul = max(ten,min(hundred,eps**meighth))
           tol = tolmul*eps
           thresh = max(tol,maxitr*q*q*unfl)
           ! test for negligible sines or cosines
           do i = 1,q
              if (theta(i) < thresh) then
                 theta(i) = zero
              else if (theta(i) > piover2 - thresh) then
                 theta(i) = piover2
              end if
           end do
           do i = 1,q - 1
              if (phi(i) < thresh) then
                 phi(i) = zero
              else if (phi(i) > piover2 - thresh) then
                 phi(i) = piover2
              end if
           end do
           ! initial deflation
           imax = q
           do while (imax > 1)
              if (phi(imax - 1) /= zero) then
                 exit
              end if
              imax = imax - 1
           end do
           imin = imax - 1
           if (imin > 1) then
              do while (phi(imin - 1) /= zero)
                 imin = imin - 1
                 if (imin <= 1) exit
              end do
           end if
           ! initialize iteration counter
           maxit = maxitr*q*q
           iter = 0
           ! begin main iteration loop
           do while (imax > 1)
              ! compute the matrix entries
              b11d(imin) = cos(theta(imin))
              b21d(imin) = -sin(theta(imin))
              do i = imin,imax - 1
                 b11e(i) = -sin(theta(i))*sin(phi(i))
                 b11d(i + 1) = cos(theta(i + 1))*cos(phi(i))
                 b12d(i) = sin(theta(i))*cos(phi(i))
                 b12e(i) = cos(theta(i + 1))*sin(phi(i))
                 b21e(i) = -cos(theta(i))*sin(phi(i))
                 b21d(i + 1) = -sin(theta(i + 1))*cos(phi(i))
                 b22d(i) = cos(theta(i))*cos(phi(i))
                 b22e(i) = -sin(theta(i + 1))*sin(phi(i))
              end do
              b12d(imax) = sin(theta(imax))
              b22d(imax) = cos(theta(imax))
              ! abort if not converging; otherwise, increment iter
              if (iter > maxit) then
                 info = 0
                 do i = 1,q
                    if (phi(i) /= zero) info = info + 1
                 end do
                 return
              end if
              iter = iter + imax - imin
              ! compute shifts
              thetamax = theta(imin)
              thetamin = theta(imin)
              do i = imin + 1,imax
                 if (theta(i) > thetamax) thetamax = theta(i)
                 if (theta(i) < thetamin) thetamin = theta(i)
              end do
              if (thetamax > piover2 - thresh) then
                 ! zero on diagonals of b11 and b22; induce deflation with a
                 ! zero shift
                 mu = zero
                 nu = one
              else if (thetamin < thresh) then
                 ! zero on diagonals of b12 and b22; induce deflation with a
                 ! zero shift
                 mu = one
                 nu = zero
              else
                 ! compute shifts for b11 and b21 and use the lesser
                 call la_slas2(b11d(imax - 1),b11e(imax - 1),b11d(imax),sigma11,dummy)

                 call la_slas2(b21d(imax - 1),b21e(imax - 1),b21d(imax),sigma21,dummy)

                 if (sigma11 <= sigma21) then
                    mu = sigma11
                    nu = sqrt(one - mu**2)
                    if (mu < thresh) then
                       mu = zero
                       nu = one
                    end if
                 else
                    nu = sigma21
                    mu = sqrt(1.0_sp - nu**2)
                    if (nu < thresh) then
                       mu = one
                       nu = zero
                    end if
                 end if
              end if
              ! rotate to produce bulges in b11 and b21
              if (mu <= nu) then
                 call la_slartgs(b11d(imin),b11e(imin),mu,rwork(iv1tcs + imin - 1),rwork( &
                           iv1tsn + imin - 1))
              else
                 call la_slartgs(b21d(imin),b21e(imin),nu,rwork(iv1tcs + imin - 1),rwork( &
                           iv1tsn + imin - 1))
              end if
              temp = rwork(iv1tcs + imin - 1)*b11d(imin) + rwork(iv1tsn + imin - 1)*b11e(imin)
              b11e(imin) = rwork(iv1tcs + imin - 1)*b11e(imin) - rwork(iv1tsn + imin - 1)*b11d(imin)

              b11d(imin) = temp
              b11bulge = rwork(iv1tsn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = rwork(iv1tcs + imin - 1)*b11d(imin + 1)
              temp = rwork(iv1tcs + imin - 1)*b21d(imin) + rwork(iv1tsn + imin - 1)*b21e(imin)
              b21e(imin) = rwork(iv1tcs + imin - 1)*b21e(imin) - rwork(iv1tsn + imin - 1)*b21d(imin)

              b21d(imin) = temp
              b21bulge = rwork(iv1tsn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = rwork(iv1tcs + imin - 1)*b21d(imin + 1)
              ! compute theta(imin)
              theta(imin) = atan2(sqrt(b21d(imin)**2 + b21bulge**2),sqrt(b11d(imin)**2 + &
                        b11bulge**2))
              ! chase the bulges in b11(imin+1,imin) and b21(imin+1,imin)
              if (b11d(imin)**2 + b11bulge**2 > thresh**2) then
                 call la_slartgp(b11bulge,b11d(imin),rwork(iu1sn + imin - 1),rwork(iu1cs + imin - &
                           1),r)
              else if (mu <= nu) then
                 call la_slartgs(b11e(imin),b11d(imin + 1),mu,rwork(iu1cs + imin - 1), &
                           rwork(iu1sn + imin - 1))
              else
                 call la_slartgs(b12d(imin),b12e(imin),nu,rwork(iu1cs + imin - 1),rwork( &
                           iu1sn + imin - 1))
              end if
              if (b21d(imin)**2 + b21bulge**2 > thresh**2) then
                 call la_slartgp(b21bulge,b21d(imin),rwork(iu2sn + imin - 1),rwork(iu2cs + imin - &
                           1),r)
              else if (nu < mu) then
                 call la_slartgs(b21e(imin),b21d(imin + 1),nu,rwork(iu2cs + imin - 1), &
                           rwork(iu2sn + imin - 1))
              else
                 call la_slartgs(b22d(imin),b22e(imin),mu,rwork(iu2cs + imin - 1),rwork(iu2sn + &
                           imin - 1))
              end if
              rwork(iu2cs + imin - 1) = -rwork(iu2cs + imin - 1)
              rwork(iu2sn + imin - 1) = -rwork(iu2sn + imin - 1)
              temp = rwork(iu1cs + imin - 1)*b11e(imin) + rwork(iu1sn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = rwork(iu1cs + imin - 1)*b11d(imin + 1) - rwork(iu1sn + imin - 1)*b11e(imin)

              b11e(imin) = temp
              if (imax > imin + 1) then
                 b11bulge = rwork(iu1sn + imin - 1)*b11e(imin + 1)
                 b11e(imin + 1) = rwork(iu1cs + imin - 1)*b11e(imin + 1)
              end if
              temp = rwork(iu1cs + imin - 1)*b12d(imin) + rwork(iu1sn + imin - 1)*b12e(imin)
              b12e(imin) = rwork(iu1cs + imin - 1)*b12e(imin) - rwork(iu1sn + imin - 1)*b12d(imin)
              b12d(imin) = temp
              b12bulge = rwork(iu1sn + imin - 1)*b12d(imin + 1)
              b12d(imin + 1) = rwork(iu1cs + imin - 1)*b12d(imin + 1)
              temp = rwork(iu2cs + imin - 1)*b21e(imin) + rwork(iu2sn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = rwork(iu2cs + imin - 1)*b21d(imin + 1) - rwork(iu2sn + imin - 1)*b21e(imin)

              b21e(imin) = temp
              if (imax > imin + 1) then
                 b21bulge = rwork(iu2sn + imin - 1)*b21e(imin + 1)
                 b21e(imin + 1) = rwork(iu2cs + imin - 1)*b21e(imin + 1)
              end if
              temp = rwork(iu2cs + imin - 1)*b22d(imin) + rwork(iu2sn + imin - 1)*b22e(imin)
              b22e(imin) = rwork(iu2cs + imin - 1)*b22e(imin) - rwork(iu2sn + imin - 1)*b22d(imin)
              b22d(imin) = temp
              b22bulge = rwork(iu2sn + imin - 1)*b22d(imin + 1)
              b22d(imin + 1) = rwork(iu2cs + imin - 1)*b22d(imin + 1)
              ! inner loop: chase bulges from b11(imin,imin+2),
              ! b12(imin,imin+1), b21(imin,imin+2), and b22(imin,imin+1) to
              ! bottom-right
              do i = imin + 1,imax - 1
                 ! compute phi(i-1)
                 x1 = sin(theta(i - 1))*b11e(i - 1) + cos(theta(i - 1))*b21e(i - 1)
                 x2 = sin(theta(i - 1))*b11bulge + cos(theta(i - 1))*b21bulge
                 y1 = sin(theta(i - 1))*b12d(i - 1) + cos(theta(i - 1))*b22d(i - 1)
                 y2 = sin(theta(i - 1))*b12bulge + cos(theta(i - 1))*b22bulge
                 phi(i - 1) = atan2(sqrt(x1**2 + x2**2),sqrt(y1**2 + y2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11e(i - 1)**2 + b11bulge**2 <= thresh**2
                 restart21 = b21e(i - 1)**2 + b21bulge**2 <= thresh**2
                 restart12 = b12d(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart22 = b22d(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i-1,i+1), b12(i-1,i),
                 ! b21(i-1,i+1), and b22(i-1,i). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart21) then
                    call la_slartgp(x2,x1,rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1),r)
                 else if (.not. restart11 .and. restart21) then
                    call la_slartgp(b11bulge,b11e(i - 1),rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1), &
                               r)
                 else if (restart11 .and. .not. restart21) then
                    call la_slartgp(b21bulge,b21e(i - 1),rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1), &
                               r)
                 else if (mu <= nu) then
                    call la_slartgs(b11d(i),b11e(i),mu,rwork(iv1tcs + i - 1),rwork(iv1tsn + i - 1) &
                               )
                 else
                    call la_slartgs(b21d(i),b21e(i),nu,rwork(iv1tcs + i - 1),rwork(iv1tsn + i - 1) &
                               )
                 end if
                 rwork(iv1tcs + i - 1) = -rwork(iv1tcs + i - 1)
                 rwork(iv1tsn + i - 1) = -rwork(iv1tsn + i - 1)
                 if (.not. restart12 .and. .not. restart22) then
                    call la_slartgp(y2,y1,rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - 1 - 1),r)

                 else if (.not. restart12 .and. restart22) then
                    call la_slartgp(b12bulge,b12d(i - 1),rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - &
                              1 - 1),r)
                 else if (restart12 .and. .not. restart22) then
                    call la_slartgp(b22bulge,b22d(i - 1),rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - &
                              1 - 1),r)
                 else if (nu < mu) then
                    call la_slartgs(b12e(i - 1),b12d(i),nu,rwork(iv2tcs + i - 1 - 1),rwork(iv2tsn + &
                              i - 1 - 1))
                 else
                    call la_slartgs(b22e(i - 1),b22d(i),mu,rwork(iv2tcs + i - 1 - 1),rwork(iv2tsn + &
                              i - 1 - 1))
                 end if
                 temp = rwork(iv1tcs + i - 1)*b11d(i) + rwork(iv1tsn + i - 1)*b11e(i)
                 b11e(i) = rwork(iv1tcs + i - 1)*b11e(i) - rwork(iv1tsn + i - 1)*b11d(i)
                 b11d(i) = temp
                 b11bulge = rwork(iv1tsn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = rwork(iv1tcs + i - 1)*b11d(i + 1)
                 temp = rwork(iv1tcs + i - 1)*b21d(i) + rwork(iv1tsn + i - 1)*b21e(i)
                 b21e(i) = rwork(iv1tcs + i - 1)*b21e(i) - rwork(iv1tsn + i - 1)*b21d(i)
                 b21d(i) = temp
                 b21bulge = rwork(iv1tsn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = rwork(iv1tcs + i - 1)*b21d(i + 1)
                 temp = rwork(iv2tcs + i - 1 - 1)*b12e(i - 1) + rwork(iv2tsn + i - 1 - 1)*b12d(i)
                 b12d(i) = rwork(iv2tcs + i - 1 - 1)*b12d(i) - rwork(iv2tsn + i - 1 - 1)*b12e(i - 1)
                 b12e(i - 1) = temp
                 b12bulge = rwork(iv2tsn + i - 1 - 1)*b12e(i)
                 b12e(i) = rwork(iv2tcs + i - 1 - 1)*b12e(i)
                 temp = rwork(iv2tcs + i - 1 - 1)*b22e(i - 1) + rwork(iv2tsn + i - 1 - 1)*b22d(i)
                 b22d(i) = rwork(iv2tcs + i - 1 - 1)*b22d(i) - rwork(iv2tsn + i - 1 - 1)*b22e(i - 1)
                 b22e(i - 1) = temp
                 b22bulge = rwork(iv2tsn + i - 1 - 1)*b22e(i)
                 b22e(i) = rwork(iv2tcs + i - 1 - 1)*b22e(i)
                 ! compute theta(i)
                 x1 = cos(phi(i - 1))*b11d(i) + sin(phi(i - 1))*b12e(i - 1)
                 x2 = cos(phi(i - 1))*b11bulge + sin(phi(i - 1))*b12bulge
                 y1 = cos(phi(i - 1))*b21d(i) + sin(phi(i - 1))*b22e(i - 1)
                 y2 = cos(phi(i - 1))*b21bulge + sin(phi(i - 1))*b22bulge
                 theta(i) = atan2(sqrt(y1**2 + y2**2),sqrt(x1**2 + x2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11d(i)**2 + b11bulge**2 <= thresh**2
                 restart12 = b12e(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart21 = b21d(i)**2 + b21bulge**2 <= thresh**2
                 restart22 = b22e(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i+1,i), b12(i+1,i-1),
                 ! b21(i+1,i), and b22(i+1,i-1). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart12) then
                    call la_slartgp(x2,x1,rwork(iu1sn + i - 1),rwork(iu1cs + i - 1),r)
                 else if (.not. restart11 .and. restart12) then
                    call la_slartgp(b11bulge,b11d(i),rwork(iu1sn + i - 1),rwork(iu1cs + i - 1),r)

                 else if (restart11 .and. .not. restart12) then
                    call la_slartgp(b12bulge,b12e(i - 1),rwork(iu1sn + i - 1),rwork(iu1cs + i - 1), &
                              r)
                 else if (mu <= nu) then
                    call la_slartgs(b11e(i),b11d(i + 1),mu,rwork(iu1cs + i - 1),rwork(iu1sn + i - 1) &
                               )
                 else
                    call la_slartgs(b12d(i),b12e(i),nu,rwork(iu1cs + i - 1),rwork(iu1sn + i - 1))

                 end if
                 if (.not. restart21 .and. .not. restart22) then
                    call la_slartgp(y2,y1,rwork(iu2sn + i - 1),rwork(iu2cs + i - 1),r)
                 else if (.not. restart21 .and. restart22) then
                    call la_slartgp(b21bulge,b21d(i),rwork(iu2sn + i - 1),rwork(iu2cs + i - 1),r)

                 else if (restart21 .and. .not. restart22) then
                    call la_slartgp(b22bulge,b22e(i - 1),rwork(iu2sn + i - 1),rwork(iu2cs + i - 1), &
                              r)
                 else if (nu < mu) then
                    call la_slartgs(b21e(i),b21e(i + 1),nu,rwork(iu2cs + i - 1),rwork(iu2sn + i - 1) &
                               )
                 else
                    call la_slartgs(b22d(i),b22e(i),mu,rwork(iu2cs + i - 1),rwork(iu2sn + i - 1))

                 end if
                 rwork(iu2cs + i - 1) = -rwork(iu2cs + i - 1)
                 rwork(iu2sn + i - 1) = -rwork(iu2sn + i - 1)
                 temp = rwork(iu1cs + i - 1)*b11e(i) + rwork(iu1sn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = rwork(iu1cs + i - 1)*b11d(i + 1) - rwork(iu1sn + i - 1)*b11e(i)
                 b11e(i) = temp
                 if (i < imax - 1) then
                    b11bulge = rwork(iu1sn + i - 1)*b11e(i + 1)
                    b11e(i + 1) = rwork(iu1cs + i - 1)*b11e(i + 1)
                 end if
                 temp = rwork(iu2cs + i - 1)*b21e(i) + rwork(iu2sn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = rwork(iu2cs + i - 1)*b21d(i + 1) - rwork(iu2sn + i - 1)*b21e(i)
                 b21e(i) = temp
                 if (i < imax - 1) then
                    b21bulge = rwork(iu2sn + i - 1)*b21e(i + 1)
                    b21e(i + 1) = rwork(iu2cs + i - 1)*b21e(i + 1)
                 end if
                 temp = rwork(iu1cs + i - 1)*b12d(i) + rwork(iu1sn + i - 1)*b12e(i)
                 b12e(i) = rwork(iu1cs + i - 1)*b12e(i) - rwork(iu1sn + i - 1)*b12d(i)
                 b12d(i) = temp
                 b12bulge = rwork(iu1sn + i - 1)*b12d(i + 1)
                 b12d(i + 1) = rwork(iu1cs + i - 1)*b12d(i + 1)
                 temp = rwork(iu2cs + i - 1)*b22d(i) + rwork(iu2sn + i - 1)*b22e(i)
                 b22e(i) = rwork(iu2cs + i - 1)*b22e(i) - rwork(iu2sn + i - 1)*b22d(i)
                 b22d(i) = temp
                 b22bulge = rwork(iu2sn + i - 1)*b22d(i + 1)
                 b22d(i + 1) = rwork(iu2cs + i - 1)*b22d(i + 1)
              end do
              ! compute phi(imax-1)
              x1 = sin(theta(imax - 1))*b11e(imax - 1) + cos(theta(imax - 1))*b21e(imax - 1)
              y1 = sin(theta(imax - 1))*b12d(imax - 1) + cos(theta(imax - 1))*b22d(imax - 1)
              y2 = sin(theta(imax - 1))*b12bulge + cos(theta(imax - 1))*b22bulge
              phi(imax - 1) = atan2(abs(x1),sqrt(y1**2 + y2**2))
              ! chase bulges from b12(imax-1,imax) and b22(imax-1,imax)
              restart12 = b12d(imax - 1)**2 + b12bulge**2 <= thresh**2
              restart22 = b22d(imax - 1)**2 + b22bulge**2 <= thresh**2
              if (.not. restart12 .and. .not. restart22) then
                 call la_slartgp(y2,y1,rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + imax - 1 - 1),r)

              else if (.not. restart12 .and. restart22) then
                 call la_slartgp(b12bulge,b12d(imax - 1),rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + &
                           imax - 1 - 1),r)
              else if (restart12 .and. .not. restart22) then
                 call la_slartgp(b22bulge,b22d(imax - 1),rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + &
                           imax - 1 - 1),r)
              else if (nu < mu) then
                 call la_slartgs(b12e(imax - 1),b12d(imax),nu,rwork(iv2tcs + imax - 1 - 1),rwork( &
                           iv2tsn + imax - 1 - 1))
              else
                 call la_slartgs(b22e(imax - 1),b22d(imax),mu,rwork(iv2tcs + imax - 1 - 1),rwork( &
                           iv2tsn + imax - 1 - 1))
              end if
              temp = rwork(iv2tcs + imax - 1 - 1)*b12e(imax - 1) + rwork(iv2tsn + imax - 1 - 1)*b12d(imax)

              b12d(imax) = rwork(iv2tcs + imax - 1 - 1)*b12d(imax) - rwork(iv2tsn + imax - 1 - 1)*b12e(imax - 1)

              b12e(imax - 1) = temp
              temp = rwork(iv2tcs + imax - 1 - 1)*b22e(imax - 1) + rwork(iv2tsn + imax - 1 - 1)*b22d(imax)

              b22d(imax) = rwork(iv2tcs + imax - 1 - 1)*b22d(imax) - rwork(iv2tsn + imax - 1 - 1)*b22e(imax - 1)

              b22e(imax - 1) = temp
              ! update singular vectors
              if (wantu1) then
                 if (colmajor) then
                    call la_clasr('R','V','F',p,imax - imin + 1,rwork(iu1cs + imin - 1),rwork( &
                              iu1sn + imin - 1),u1(1,imin),ldu1)
                 else
                    call la_clasr('L','V','F',imax - imin + 1,p,rwork(iu1cs + imin - 1),rwork( &
                              iu1sn + imin - 1),u1(imin,1),ldu1)
                 end if
              end if
              if (wantu2) then
                 if (colmajor) then
                    call la_clasr('R','V','F',m - p,imax - imin + 1,rwork(iu2cs + imin - 1),rwork( &
                              iu2sn + imin - 1),u2(1,imin),ldu2)
                 else
                    call la_clasr('L','V','F',imax - imin + 1,m - p,rwork(iu2cs + imin - 1),rwork( &
                              iu2sn + imin - 1),u2(imin,1),ldu2)
                 end if
              end if
              if (wantv1t) then
                 if (colmajor) then
                    call la_clasr('L','V','F',imax - imin + 1,q,rwork(iv1tcs + imin - 1),rwork( &
                              iv1tsn + imin - 1),v1t(imin,1),ldv1t)
                 else
                    call la_clasr('R','V','F',q,imax - imin + 1,rwork(iv1tcs + imin - 1),rwork( &
                              iv1tsn + imin - 1),v1t(1,imin),ldv1t)
                 end if
              end if
              if (wantv2t) then
                 if (colmajor) then
                    call la_clasr('L','V','F',imax - imin + 1,m - q,rwork(iv2tcs + imin - 1), &
                              rwork(iv2tsn + imin - 1),v2t(imin,1),ldv2t)
                 else
                    call la_clasr('R','V','F',m - q,imax - imin + 1,rwork(iv2tcs + imin - 1), &
                              rwork(iv2tsn + imin - 1),v2t(1,imin),ldv2t)
                 end if
              end if
              ! fix signs on b11(imax-1,imax) and b21(imax-1,imax)
              if (b11e(imax - 1) + b21e(imax - 1) > 0) then
                 b11d(imax) = -b11d(imax)
                 b21d(imax) = -b21d(imax)
                 if (wantv1t) then
                    if (colmajor) then
                       call la_cscal(q,cnegone,v1t(imax,1),ldv1t)
                    else
                       call la_cscal(q,cnegone,v1t(1,imax),1)
                    end if
                 end if
              end if
              ! compute theta(imax)
              x1 = cos(phi(imax - 1))*b11d(imax) + sin(phi(imax - 1))*b12e(imax - 1)
              y1 = cos(phi(imax - 1))*b21d(imax) + sin(phi(imax - 1))*b22e(imax - 1)
              theta(imax) = atan2(abs(y1),abs(x1))
              ! fix signs on b11(imax,imax), b12(imax,imax-1), b21(imax,imax),
              ! and b22(imax,imax-1)
              if (b11d(imax) + b12e(imax - 1) < 0) then
                 b12d(imax) = -b12d(imax)
                 if (wantu1) then
                    if (colmajor) then
                       call la_cscal(p,cnegone,u1(1,imax),1)
                    else
                       call la_cscal(p,cnegone,u1(imax,1),ldu1)
                    end if
                 end if
              end if
              if (b21d(imax) + b22e(imax - 1) > 0) then
                 b22d(imax) = -b22d(imax)
                 if (wantu2) then
                    if (colmajor) then
                       call la_cscal(m - p,cnegone,u2(1,imax),1)
                    else
                       call la_cscal(m - p,cnegone,u2(imax,1),ldu2)
                    end if
                 end if
              end if
              ! fix signs on b12(imax,imax) and b22(imax,imax)
              if (b12d(imax) + b22d(imax) < 0) then
                 if (wantv2t) then
                    if (colmajor) then
                       call la_cscal(m - q,cnegone,v2t(imax,1),ldv2t)
                    else
                       call la_cscal(m - q,cnegone,v2t(1,imax),1)
                    end if
                 end if
              end if
              ! test for negligible sines or cosines
              do i = imin,imax
                 if (theta(i) < thresh) then
                    theta(i) = zero
                 else if (theta(i) > piover2 - thresh) then
                    theta(i) = piover2
                 end if
              end do
              do i = imin,imax - 1
                 if (phi(i) < thresh) then
                    phi(i) = zero
                 else if (phi(i) > piover2 - thresh) then
                    phi(i) = piover2
                 end if
              end do
              ! deflate
              if (imax > 1) then
                 do while (phi(imax - 1) == zero)
                    imax = imax - 1
                    if (imax <= 1) exit
                 end do
              end if
              if (imin > imax - 1) imin = imax - 1
              if (imin > 1) then
                 do while (phi(imin - 1) /= zero)
                     imin = imin - 1
                     if (imin <= 1) exit
                 end do
              end if
              ! repeat main iteration loop
           end do
           ! postprocessing: order theta from least to greatest
           do i = 1,q
              mini = i
              thetamin = theta(i)
              do j = i + 1,q
                 if (theta(j) < thetamin) then
                    mini = j
                    thetamin = theta(j)
                 end if
              end do
              if (mini /= i) then
                 theta(mini) = theta(i)
                 theta(i) = thetamin
                 if (colmajor) then
                    if (wantu1) call la_cswap(p,u1(1,i),1,u1(1,mini),1)
                    if (wantu2) call la_cswap(m - p,u2(1,i),1,u2(1,mini),1)
                    if (wantv1t) call la_cswap(q,v1t(i,1),ldv1t,v1t(mini,1),ldv1t)

                    if (wantv2t) call la_cswap(m - q,v2t(i,1),ldv2t,v2t(mini,1),ldv2t)

                 else
                    if (wantu1) call la_cswap(p,u1(i,1),ldu1,u1(mini,1),ldu1)
                    if (wantu2) call la_cswap(m - p,u2(i,1),ldu2,u2(mini,1),ldu2)
                    if (wantv1t) call la_cswap(q,v1t(1,i),1,v1t(1,mini),1)
                    if (wantv2t) call la_cswap(m - q,v2t(1,i),1,v2t(1,mini),1)
                 end if
              end if
           end do
           return
     end subroutine la_cbbcsd
     !> ZBBCSD: computes the CS decomposition of a unitary matrix in
     !> bidiagonal-block form,
     !> [ B11 | B12 0  0 ]
     !> [  0  |  0 -I  0 ]
     !> X = [----------------]
     !> [ B21 | B22 0  0 ]
     !> [  0  |  0  0  I ]
     !> [  C | -S  0  0 ]
     !> [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**H
     !> = [---------] [---------------] [---------]   .
     !> [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
     !> [  0 |  0  0  I ]
     !> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
     !> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
     !> transposed and/or permuted. This can be done in constant time using
     !> the TRANS and SIGNS options. See ZUNCSD for details.)
     !> The bidiagonal matrices B11, B12, B21, and B22 are represented
     !> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
     !> The unitary matrices U1, U2, V1T, and V2T are input/output.
     !> The input matrices are pre- or post-multiplied by the appropriate
     !> singular vector matrices.

     pure subroutine la_zbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,phi,u1, &
     ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,b11d,b11e,b12d,b12e,b21d,b21e,b22d,b22e,rwork, &
               lrwork,info)
        use la_constants_dp,only:zero,one,ten,cnegone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,lrwork,m,p,q
           ! Array Arguments
           real(dp),intent(out) :: b11d(*),b11e(*),b12d(*),b12e(*),b21d(*),b21e(*),b22d(*), &
                      b22e(*),rwork(*)
           real(dp),intent(inout) :: phi(*),theta(*)
           complex(dp),intent(inout) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*)

        ! ===================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 6
           real(dp),parameter :: hundred = 100.0_dp
           real(dp),parameter :: meighth = -0.125_dp
           real(dp),parameter :: piover2 = 1.57079632679489661923132169163975144210_dp

           ! Local Scalars
           logical(lk) :: colmajor,lquery,restart11,restart12,restart21,restart22,wantu1, &
                     wantu2,wantv1t,wantv2t
           integer(ilp) :: i,imin,imax,iter,iu1cs,iu1sn,iu2cs,iu2sn,iv1tcs,iv1tsn, &
                     iv2tcs,iv2tsn,j,lrworkmin,lrworkopt,maxit,mini
           real(dp) :: b11bulge,b12bulge,b21bulge,b22bulge,dummy,eps,mu,nu,r,sigma11, &
                     sigma21,temp,thetamax,thetamin,thresh,tol,tolmul,unfl,x1,x2,y1,y2
           ! Intrinsic Functions
           intrinsic :: abs,atan2,cos,max,min,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lrwork == -1
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           if (m < 0) then
              info = -6
           else if (p < 0 .or. p > m) then
              info = -7
           else if (q < 0 .or. q > m) then
              info = -8
           else if (q > p .or. q > m - p .or. q > m - q) then
              info = -8
           else if (wantu1 .and. ldu1 < p) then
              info = -12
           else if (wantu2 .and. ldu2 < m - p) then
              info = -14
           else if (wantv1t .and. ldv1t < q) then
              info = -16
           else if (wantv2t .and. ldv2t < m - q) then
              info = -18
           end if
           ! quick return if q = 0
           if (info == 0 .and. q == 0) then
              lrworkmin = 1
              rwork(1) = lrworkmin
              return
           end if
           ! compute workspace
           if (info == 0) then
              iu1cs = 1
              iu1sn = iu1cs + q
              iu2cs = iu1sn + q
              iu2sn = iu2cs + q
              iv1tcs = iu2sn + q
              iv1tsn = iv1tcs + q
              iv2tcs = iv1tsn + q
              iv2tsn = iv2tcs + q
              lrworkopt = iv2tsn + q - 1
              lrworkmin = lrworkopt
              rwork(1) = lrworkopt
              if (lrwork < lrworkmin .and. .not. lquery) then
                 info = -28
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZBBCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! get machine constants
           eps = la_dlamch('EPSILON')
           unfl = la_dlamch('SAFE MINIMUM')
           tolmul = max(ten,min(hundred,eps**meighth))
           tol = tolmul*eps
           thresh = max(tol,maxitr*q*q*unfl)
           ! test for negligible sines or cosines
           do i = 1,q
              if (theta(i) < thresh) then
                 theta(i) = zero
              else if (theta(i) > piover2 - thresh) then
                 theta(i) = piover2
              end if
           end do
           do i = 1,q - 1
              if (phi(i) < thresh) then
                 phi(i) = zero
              else if (phi(i) > piover2 - thresh) then
                 phi(i) = piover2
              end if
           end do
           ! initial deflation
           imax = q
           do while (imax > 1)
              if (phi(imax - 1) /= zero) then
                 exit
              end if
              imax = imax - 1
           end do
           imin = imax - 1
           if (imin > 1) then
              do while (phi(imin - 1) /= zero)
                 imin = imin - 1
                 if (imin <= 1) exit
              end do
           end if
           ! initialize iteration counter
           maxit = maxitr*q*q
           iter = 0
           ! begin main iteration loop
           do while (imax > 1)
              ! compute the matrix entries
              b11d(imin) = cos(theta(imin))
              b21d(imin) = -sin(theta(imin))
              do i = imin,imax - 1
                 b11e(i) = -sin(theta(i))*sin(phi(i))
                 b11d(i + 1) = cos(theta(i + 1))*cos(phi(i))
                 b12d(i) = sin(theta(i))*cos(phi(i))
                 b12e(i) = cos(theta(i + 1))*sin(phi(i))
                 b21e(i) = -cos(theta(i))*sin(phi(i))
                 b21d(i + 1) = -sin(theta(i + 1))*cos(phi(i))
                 b22d(i) = cos(theta(i))*cos(phi(i))
                 b22e(i) = -sin(theta(i + 1))*sin(phi(i))
              end do
              b12d(imax) = sin(theta(imax))
              b22d(imax) = cos(theta(imax))
              ! abort if not converging; otherwise, increment iter
              if (iter > maxit) then
                 info = 0
                 do i = 1,q
                    if (phi(i) /= zero) info = info + 1
                 end do
                 return
              end if
              iter = iter + imax - imin
              ! compute shifts
              thetamax = theta(imin)
              thetamin = theta(imin)
              do i = imin + 1,imax
                 if (theta(i) > thetamax) thetamax = theta(i)
                 if (theta(i) < thetamin) thetamin = theta(i)
              end do
              if (thetamax > piover2 - thresh) then
                 ! zero on diagonals of b11 and b22; induce deflation with a
                 ! zero shift
                 mu = zero
                 nu = one
              else if (thetamin < thresh) then
                 ! zero on diagonals of b12 and b22; induce deflation with a
                 ! zero shift
                 mu = one
                 nu = zero
              else
                 ! compute shifts for b11 and b21 and use the lesser
                 call la_dlas2(b11d(imax - 1),b11e(imax - 1),b11d(imax),sigma11,dummy)

                 call la_dlas2(b21d(imax - 1),b21e(imax - 1),b21d(imax),sigma21,dummy)

                 if (sigma11 <= sigma21) then
                    mu = sigma11
                    nu = sqrt(one - mu**2)
                    if (mu < thresh) then
                       mu = zero
                       nu = one
                    end if
                 else
                    nu = sigma21
                    mu = sqrt(1.0_dp - nu**2)
                    if (nu < thresh) then
                       mu = one
                       nu = zero
                    end if
                 end if
              end if
              ! rotate to produce bulges in b11 and b21
              if (mu <= nu) then
                 call la_dlartgs(b11d(imin),b11e(imin),mu,rwork(iv1tcs + imin - 1),rwork( &
                           iv1tsn + imin - 1))
              else
                 call la_dlartgs(b21d(imin),b21e(imin),nu,rwork(iv1tcs + imin - 1),rwork( &
                           iv1tsn + imin - 1))
              end if
              temp = rwork(iv1tcs + imin - 1)*b11d(imin) + rwork(iv1tsn + imin - 1)*b11e(imin)
              b11e(imin) = rwork(iv1tcs + imin - 1)*b11e(imin) - rwork(iv1tsn + imin - 1)*b11d(imin)

              b11d(imin) = temp
              b11bulge = rwork(iv1tsn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = rwork(iv1tcs + imin - 1)*b11d(imin + 1)
              temp = rwork(iv1tcs + imin - 1)*b21d(imin) + rwork(iv1tsn + imin - 1)*b21e(imin)
              b21e(imin) = rwork(iv1tcs + imin - 1)*b21e(imin) - rwork(iv1tsn + imin - 1)*b21d(imin)

              b21d(imin) = temp
              b21bulge = rwork(iv1tsn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = rwork(iv1tcs + imin - 1)*b21d(imin + 1)
              ! compute theta(imin)
              theta(imin) = atan2(sqrt(b21d(imin)**2 + b21bulge**2),sqrt(b11d(imin)**2 + &
                        b11bulge**2))
              ! chase the bulges in b11(imin+1,imin) and b21(imin+1,imin)
              if (b11d(imin)**2 + b11bulge**2 > thresh**2) then
                 call la_dlartgp(b11bulge,b11d(imin),rwork(iu1sn + imin - 1),rwork(iu1cs + imin - &
                           1),r)
              else if (mu <= nu) then
                 call la_dlartgs(b11e(imin),b11d(imin + 1),mu,rwork(iu1cs + imin - 1), &
                           rwork(iu1sn + imin - 1))
              else
                 call la_dlartgs(b12d(imin),b12e(imin),nu,rwork(iu1cs + imin - 1),rwork( &
                           iu1sn + imin - 1))
              end if
              if (b21d(imin)**2 + b21bulge**2 > thresh**2) then
                 call la_dlartgp(b21bulge,b21d(imin),rwork(iu2sn + imin - 1),rwork(iu2cs + imin - &
                           1),r)
              else if (nu < mu) then
                 call la_dlartgs(b21e(imin),b21d(imin + 1),nu,rwork(iu2cs + imin - 1), &
                           rwork(iu2sn + imin - 1))
              else
                 call la_dlartgs(b22d(imin),b22e(imin),mu,rwork(iu2cs + imin - 1),rwork(iu2sn + &
                           imin - 1))
              end if
              rwork(iu2cs + imin - 1) = -rwork(iu2cs + imin - 1)
              rwork(iu2sn + imin - 1) = -rwork(iu2sn + imin - 1)
              temp = rwork(iu1cs + imin - 1)*b11e(imin) + rwork(iu1sn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = rwork(iu1cs + imin - 1)*b11d(imin + 1) - rwork(iu1sn + imin - 1)*b11e(imin)

              b11e(imin) = temp
              if (imax > imin + 1) then
                 b11bulge = rwork(iu1sn + imin - 1)*b11e(imin + 1)
                 b11e(imin + 1) = rwork(iu1cs + imin - 1)*b11e(imin + 1)
              end if
              temp = rwork(iu1cs + imin - 1)*b12d(imin) + rwork(iu1sn + imin - 1)*b12e(imin)
              b12e(imin) = rwork(iu1cs + imin - 1)*b12e(imin) - rwork(iu1sn + imin - 1)*b12d(imin)
              b12d(imin) = temp
              b12bulge = rwork(iu1sn + imin - 1)*b12d(imin + 1)
              b12d(imin + 1) = rwork(iu1cs + imin - 1)*b12d(imin + 1)
              temp = rwork(iu2cs + imin - 1)*b21e(imin) + rwork(iu2sn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = rwork(iu2cs + imin - 1)*b21d(imin + 1) - rwork(iu2sn + imin - 1)*b21e(imin)

              b21e(imin) = temp
              if (imax > imin + 1) then
                 b21bulge = rwork(iu2sn + imin - 1)*b21e(imin + 1)
                 b21e(imin + 1) = rwork(iu2cs + imin - 1)*b21e(imin + 1)
              end if
              temp = rwork(iu2cs + imin - 1)*b22d(imin) + rwork(iu2sn + imin - 1)*b22e(imin)
              b22e(imin) = rwork(iu2cs + imin - 1)*b22e(imin) - rwork(iu2sn + imin - 1)*b22d(imin)
              b22d(imin) = temp
              b22bulge = rwork(iu2sn + imin - 1)*b22d(imin + 1)
              b22d(imin + 1) = rwork(iu2cs + imin - 1)*b22d(imin + 1)
              ! inner loop: chase bulges from b11(imin,imin+2),
              ! b12(imin,imin+1), b21(imin,imin+2), and b22(imin,imin+1) to
              ! bottom-right
              do i = imin + 1,imax - 1
                 ! compute phi(i-1)
                 x1 = sin(theta(i - 1))*b11e(i - 1) + cos(theta(i - 1))*b21e(i - 1)
                 x2 = sin(theta(i - 1))*b11bulge + cos(theta(i - 1))*b21bulge
                 y1 = sin(theta(i - 1))*b12d(i - 1) + cos(theta(i - 1))*b22d(i - 1)
                 y2 = sin(theta(i - 1))*b12bulge + cos(theta(i - 1))*b22bulge
                 phi(i - 1) = atan2(sqrt(x1**2 + x2**2),sqrt(y1**2 + y2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11e(i - 1)**2 + b11bulge**2 <= thresh**2
                 restart21 = b21e(i - 1)**2 + b21bulge**2 <= thresh**2
                 restart12 = b12d(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart22 = b22d(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i-1,i+1), b12(i-1,i),
                 ! b21(i-1,i+1), and b22(i-1,i). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart21) then
                    call la_dlartgp(x2,x1,rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1),r)
                 else if (.not. restart11 .and. restart21) then
                    call la_dlartgp(b11bulge,b11e(i - 1),rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1), &
                               r)
                 else if (restart11 .and. .not. restart21) then
                    call la_dlartgp(b21bulge,b21e(i - 1),rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1), &
                               r)
                 else if (mu <= nu) then
                    call la_dlartgs(b11d(i),b11e(i),mu,rwork(iv1tcs + i - 1),rwork(iv1tsn + i - 1) &
                               )
                 else
                    call la_dlartgs(b21d(i),b21e(i),nu,rwork(iv1tcs + i - 1),rwork(iv1tsn + i - 1) &
                               )
                 end if
                 rwork(iv1tcs + i - 1) = -rwork(iv1tcs + i - 1)
                 rwork(iv1tsn + i - 1) = -rwork(iv1tsn + i - 1)
                 if (.not. restart12 .and. .not. restart22) then
                    call la_dlartgp(y2,y1,rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - 1 - 1),r)

                 else if (.not. restart12 .and. restart22) then
                    call la_dlartgp(b12bulge,b12d(i - 1),rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - &
                              1 - 1),r)
                 else if (restart12 .and. .not. restart22) then
                    call la_dlartgp(b22bulge,b22d(i - 1),rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - &
                              1 - 1),r)
                 else if (nu < mu) then
                    call la_dlartgs(b12e(i - 1),b12d(i),nu,rwork(iv2tcs + i - 1 - 1),rwork(iv2tsn + &
                              i - 1 - 1))
                 else
                    call la_dlartgs(b22e(i - 1),b22d(i),mu,rwork(iv2tcs + i - 1 - 1),rwork(iv2tsn + &
                              i - 1 - 1))
                 end if
                 temp = rwork(iv1tcs + i - 1)*b11d(i) + rwork(iv1tsn + i - 1)*b11e(i)
                 b11e(i) = rwork(iv1tcs + i - 1)*b11e(i) - rwork(iv1tsn + i - 1)*b11d(i)
                 b11d(i) = temp
                 b11bulge = rwork(iv1tsn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = rwork(iv1tcs + i - 1)*b11d(i + 1)
                 temp = rwork(iv1tcs + i - 1)*b21d(i) + rwork(iv1tsn + i - 1)*b21e(i)
                 b21e(i) = rwork(iv1tcs + i - 1)*b21e(i) - rwork(iv1tsn + i - 1)*b21d(i)
                 b21d(i) = temp
                 b21bulge = rwork(iv1tsn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = rwork(iv1tcs + i - 1)*b21d(i + 1)
                 temp = rwork(iv2tcs + i - 1 - 1)*b12e(i - 1) + rwork(iv2tsn + i - 1 - 1)*b12d(i)
                 b12d(i) = rwork(iv2tcs + i - 1 - 1)*b12d(i) - rwork(iv2tsn + i - 1 - 1)*b12e(i - 1)
                 b12e(i - 1) = temp
                 b12bulge = rwork(iv2tsn + i - 1 - 1)*b12e(i)
                 b12e(i) = rwork(iv2tcs + i - 1 - 1)*b12e(i)
                 temp = rwork(iv2tcs + i - 1 - 1)*b22e(i - 1) + rwork(iv2tsn + i - 1 - 1)*b22d(i)
                 b22d(i) = rwork(iv2tcs + i - 1 - 1)*b22d(i) - rwork(iv2tsn + i - 1 - 1)*b22e(i - 1)
                 b22e(i - 1) = temp
                 b22bulge = rwork(iv2tsn + i - 1 - 1)*b22e(i)
                 b22e(i) = rwork(iv2tcs + i - 1 - 1)*b22e(i)
                 ! compute theta(i)
                 x1 = cos(phi(i - 1))*b11d(i) + sin(phi(i - 1))*b12e(i - 1)
                 x2 = cos(phi(i - 1))*b11bulge + sin(phi(i - 1))*b12bulge
                 y1 = cos(phi(i - 1))*b21d(i) + sin(phi(i - 1))*b22e(i - 1)
                 y2 = cos(phi(i - 1))*b21bulge + sin(phi(i - 1))*b22bulge
                 theta(i) = atan2(sqrt(y1**2 + y2**2),sqrt(x1**2 + x2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11d(i)**2 + b11bulge**2 <= thresh**2
                 restart12 = b12e(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart21 = b21d(i)**2 + b21bulge**2 <= thresh**2
                 restart22 = b22e(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i+1,i), b12(i+1,i-1),
                 ! b21(i+1,i), and b22(i+1,i-1). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart12) then
                    call la_dlartgp(x2,x1,rwork(iu1sn + i - 1),rwork(iu1cs + i - 1),r)
                 else if (.not. restart11 .and. restart12) then
                    call la_dlartgp(b11bulge,b11d(i),rwork(iu1sn + i - 1),rwork(iu1cs + i - 1),r)

                 else if (restart11 .and. .not. restart12) then
                    call la_dlartgp(b12bulge,b12e(i - 1),rwork(iu1sn + i - 1),rwork(iu1cs + i - 1), &
                              r)
                 else if (mu <= nu) then
                    call la_dlartgs(b11e(i),b11d(i + 1),mu,rwork(iu1cs + i - 1),rwork(iu1sn + i - 1) &
                               )
                 else
                    call la_dlartgs(b12d(i),b12e(i),nu,rwork(iu1cs + i - 1),rwork(iu1sn + i - 1))

                 end if
                 if (.not. restart21 .and. .not. restart22) then
                    call la_dlartgp(y2,y1,rwork(iu2sn + i - 1),rwork(iu2cs + i - 1),r)
                 else if (.not. restart21 .and. restart22) then
                    call la_dlartgp(b21bulge,b21d(i),rwork(iu2sn + i - 1),rwork(iu2cs + i - 1),r)

                 else if (restart21 .and. .not. restart22) then
                    call la_dlartgp(b22bulge,b22e(i - 1),rwork(iu2sn + i - 1),rwork(iu2cs + i - 1), &
                              r)
                 else if (nu < mu) then
                    call la_dlartgs(b21e(i),b21e(i + 1),nu,rwork(iu2cs + i - 1),rwork(iu2sn + i - 1) &
                               )
                 else
                    call la_dlartgs(b22d(i),b22e(i),mu,rwork(iu2cs + i - 1),rwork(iu2sn + i - 1))

                 end if
                 rwork(iu2cs + i - 1) = -rwork(iu2cs + i - 1)
                 rwork(iu2sn + i - 1) = -rwork(iu2sn + i - 1)
                 temp = rwork(iu1cs + i - 1)*b11e(i) + rwork(iu1sn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = rwork(iu1cs + i - 1)*b11d(i + 1) - rwork(iu1sn + i - 1)*b11e(i)
                 b11e(i) = temp
                 if (i < imax - 1) then
                    b11bulge = rwork(iu1sn + i - 1)*b11e(i + 1)
                    b11e(i + 1) = rwork(iu1cs + i - 1)*b11e(i + 1)
                 end if
                 temp = rwork(iu2cs + i - 1)*b21e(i) + rwork(iu2sn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = rwork(iu2cs + i - 1)*b21d(i + 1) - rwork(iu2sn + i - 1)*b21e(i)
                 b21e(i) = temp
                 if (i < imax - 1) then
                    b21bulge = rwork(iu2sn + i - 1)*b21e(i + 1)
                    b21e(i + 1) = rwork(iu2cs + i - 1)*b21e(i + 1)
                 end if
                 temp = rwork(iu1cs + i - 1)*b12d(i) + rwork(iu1sn + i - 1)*b12e(i)
                 b12e(i) = rwork(iu1cs + i - 1)*b12e(i) - rwork(iu1sn + i - 1)*b12d(i)
                 b12d(i) = temp
                 b12bulge = rwork(iu1sn + i - 1)*b12d(i + 1)
                 b12d(i + 1) = rwork(iu1cs + i - 1)*b12d(i + 1)
                 temp = rwork(iu2cs + i - 1)*b22d(i) + rwork(iu2sn + i - 1)*b22e(i)
                 b22e(i) = rwork(iu2cs + i - 1)*b22e(i) - rwork(iu2sn + i - 1)*b22d(i)
                 b22d(i) = temp
                 b22bulge = rwork(iu2sn + i - 1)*b22d(i + 1)
                 b22d(i + 1) = rwork(iu2cs + i - 1)*b22d(i + 1)
              end do
              ! compute phi(imax-1)
              x1 = sin(theta(imax - 1))*b11e(imax - 1) + cos(theta(imax - 1))*b21e(imax - 1)
              y1 = sin(theta(imax - 1))*b12d(imax - 1) + cos(theta(imax - 1))*b22d(imax - 1)
              y2 = sin(theta(imax - 1))*b12bulge + cos(theta(imax - 1))*b22bulge
              phi(imax - 1) = atan2(abs(x1),sqrt(y1**2 + y2**2))
              ! chase bulges from b12(imax-1,imax) and b22(imax-1,imax)
              restart12 = b12d(imax - 1)**2 + b12bulge**2 <= thresh**2
              restart22 = b22d(imax - 1)**2 + b22bulge**2 <= thresh**2
              if (.not. restart12 .and. .not. restart22) then
                 call la_dlartgp(y2,y1,rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + imax - 1 - 1),r)

              else if (.not. restart12 .and. restart22) then
                 call la_dlartgp(b12bulge,b12d(imax - 1),rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + &
                           imax - 1 - 1),r)
              else if (restart12 .and. .not. restart22) then
                 call la_dlartgp(b22bulge,b22d(imax - 1),rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + &
                           imax - 1 - 1),r)
              else if (nu < mu) then
                 call la_dlartgs(b12e(imax - 1),b12d(imax),nu,rwork(iv2tcs + imax - 1 - 1),rwork( &
                           iv2tsn + imax - 1 - 1))
              else
                 call la_dlartgs(b22e(imax - 1),b22d(imax),mu,rwork(iv2tcs + imax - 1 - 1),rwork( &
                           iv2tsn + imax - 1 - 1))
              end if
              temp = rwork(iv2tcs + imax - 1 - 1)*b12e(imax - 1) + rwork(iv2tsn + imax - 1 - 1)*b12d(imax)

              b12d(imax) = rwork(iv2tcs + imax - 1 - 1)*b12d(imax) - rwork(iv2tsn + imax - 1 - 1)*b12e(imax - 1)

              b12e(imax - 1) = temp
              temp = rwork(iv2tcs + imax - 1 - 1)*b22e(imax - 1) + rwork(iv2tsn + imax - 1 - 1)*b22d(imax)

              b22d(imax) = rwork(iv2tcs + imax - 1 - 1)*b22d(imax) - rwork(iv2tsn + imax - 1 - 1)*b22e(imax - 1)

              b22e(imax - 1) = temp
              ! update singular vectors
              if (wantu1) then
                 if (colmajor) then
                    call la_zlasr('R','V','F',p,imax - imin + 1,rwork(iu1cs + imin - 1),rwork( &
                              iu1sn + imin - 1),u1(1,imin),ldu1)
                 else
                    call la_zlasr('L','V','F',imax - imin + 1,p,rwork(iu1cs + imin - 1),rwork( &
                              iu1sn + imin - 1),u1(imin,1),ldu1)
                 end if
              end if
              if (wantu2) then
                 if (colmajor) then
                    call la_zlasr('R','V','F',m - p,imax - imin + 1,rwork(iu2cs + imin - 1),rwork( &
                              iu2sn + imin - 1),u2(1,imin),ldu2)
                 else
                    call la_zlasr('L','V','F',imax - imin + 1,m - p,rwork(iu2cs + imin - 1),rwork( &
                              iu2sn + imin - 1),u2(imin,1),ldu2)
                 end if
              end if
              if (wantv1t) then
                 if (colmajor) then
                    call la_zlasr('L','V','F',imax - imin + 1,q,rwork(iv1tcs + imin - 1),rwork( &
                              iv1tsn + imin - 1),v1t(imin,1),ldv1t)
                 else
                    call la_zlasr('R','V','F',q,imax - imin + 1,rwork(iv1tcs + imin - 1),rwork( &
                              iv1tsn + imin - 1),v1t(1,imin),ldv1t)
                 end if
              end if
              if (wantv2t) then
                 if (colmajor) then
                    call la_zlasr('L','V','F',imax - imin + 1,m - q,rwork(iv2tcs + imin - 1), &
                              rwork(iv2tsn + imin - 1),v2t(imin,1),ldv2t)
                 else
                    call la_zlasr('R','V','F',m - q,imax - imin + 1,rwork(iv2tcs + imin - 1), &
                              rwork(iv2tsn + imin - 1),v2t(1,imin),ldv2t)
                 end if
              end if
              ! fix signs on b11(imax-1,imax) and b21(imax-1,imax)
              if (b11e(imax - 1) + b21e(imax - 1) > 0) then
                 b11d(imax) = -b11d(imax)
                 b21d(imax) = -b21d(imax)
                 if (wantv1t) then
                    if (colmajor) then
                       call la_zscal(q,cnegone,v1t(imax,1),ldv1t)
                    else
                       call la_zscal(q,cnegone,v1t(1,imax),1)
                    end if
                 end if
              end if
              ! compute theta(imax)
              x1 = cos(phi(imax - 1))*b11d(imax) + sin(phi(imax - 1))*b12e(imax - 1)
              y1 = cos(phi(imax - 1))*b21d(imax) + sin(phi(imax - 1))*b22e(imax - 1)
              theta(imax) = atan2(abs(y1),abs(x1))
              ! fix signs on b11(imax,imax), b12(imax,imax-1), b21(imax,imax),
              ! and b22(imax,imax-1)
              if (b11d(imax) + b12e(imax - 1) < 0) then
                 b12d(imax) = -b12d(imax)
                 if (wantu1) then
                    if (colmajor) then
                       call la_zscal(p,cnegone,u1(1,imax),1)
                    else
                       call la_zscal(p,cnegone,u1(imax,1),ldu1)
                    end if
                 end if
              end if
              if (b21d(imax) + b22e(imax - 1) > 0) then
                 b22d(imax) = -b22d(imax)
                 if (wantu2) then
                    if (colmajor) then
                       call la_zscal(m - p,cnegone,u2(1,imax),1)
                    else
                       call la_zscal(m - p,cnegone,u2(imax,1),ldu2)
                    end if
                 end if
              end if
              ! fix signs on b12(imax,imax) and b22(imax,imax)
              if (b12d(imax) + b22d(imax) < 0) then
                 if (wantv2t) then
                    if (colmajor) then
                       call la_zscal(m - q,cnegone,v2t(imax,1),ldv2t)
                    else
                       call la_zscal(m - q,cnegone,v2t(1,imax),1)
                    end if
                 end if
              end if
              ! test for negligible sines or cosines
              do i = imin,imax
                 if (theta(i) < thresh) then
                    theta(i) = zero
                 else if (theta(i) > piover2 - thresh) then
                    theta(i) = piover2
                 end if
              end do
              do i = imin,imax - 1
                 if (phi(i) < thresh) then
                    phi(i) = zero
                 else if (phi(i) > piover2 - thresh) then
                    phi(i) = piover2
                 end if
              end do
              ! deflate
              if (imax > 1) then
                 do while (phi(imax - 1) == zero)
                    imax = imax - 1
                    if (imax <= 1) exit
                 end do
              end if
              if (imin > imax - 1) imin = imax - 1
              if (imin > 1) then
                 do while (phi(imin - 1) /= zero)
                     imin = imin - 1
                     if (imin <= 1) exit
                 end do
              end if
              ! repeat main iteration loop
           end do
           ! postprocessing: order theta from least to greatest
           do i = 1,q
              mini = i
              thetamin = theta(i)
              do j = i + 1,q
                 if (theta(j) < thetamin) then
                    mini = j
                    thetamin = theta(j)
                 end if
              end do
              if (mini /= i) then
                 theta(mini) = theta(i)
                 theta(i) = thetamin
                 if (colmajor) then
                    if (wantu1) call la_zswap(p,u1(1,i),1,u1(1,mini),1)
                    if (wantu2) call la_zswap(m - p,u2(1,i),1,u2(1,mini),1)
                    if (wantv1t) call la_zswap(q,v1t(i,1),ldv1t,v1t(mini,1),ldv1t)

                    if (wantv2t) call la_zswap(m - q,v2t(i,1),ldv2t,v2t(mini,1),ldv2t)

                 else
                    if (wantu1) call la_zswap(p,u1(i,1),ldu1,u1(mini,1),ldu1)
                    if (wantu2) call la_zswap(m - p,u2(i,1),ldu2,u2(mini,1),ldu2)
                    if (wantv1t) call la_zswap(q,v1t(1,i),1,v1t(1,mini),1)
                    if (wantv2t) call la_zswap(m - q,v2t(1,i),1,v2t(1,mini),1)
                 end if
              end if
           end do
           return
     end subroutine la_zbbcsd
     !> WBBCSD: computes the CS decomposition of a unitary matrix in
     !> bidiagonal-block form,
     !> [ B11 | B12 0  0 ]
     !> [  0  |  0 -I  0 ]
     !> X = [----------------]
     !> [ B21 | B22 0  0 ]
     !> [  0  |  0  0  I ]
     !> [  C | -S  0  0 ]
     !> [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**H
     !> = [---------] [---------------] [---------]   .
     !> [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
     !> [  0 |  0  0  I ]
     !> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
     !> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
     !> transposed and/or permuted. This can be done in constant time using
     !> the TRANS and SIGNS options. See WUNCSD for details.)
     !> The bidiagonal matrices B11, B12, B21, and B22 are represented
     !> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
     !> The unitary matrices U1, U2, V1T, and V2T are input/output.
     !> The input matrices are pre- or post-multiplied by the appropriate
     !> singular vector matrices.

     pure subroutine la_wbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,phi,u1, &
     ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,b11d,b11e,b12d,b12e,b21d,b21e,b22d,b22e,rwork, &
               lrwork,info)
        use la_constants_qp,only:zero,one,ten,cnegone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,lrwork,m,p,q
           ! Array Arguments
           real(qp),intent(out) :: b11d(*),b11e(*),b12d(*),b12e(*),b21d(*),b21e(*),b22d(*), &
                      b22e(*),rwork(*)
           real(qp),intent(inout) :: phi(*),theta(*)
           complex(qp),intent(inout) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*)

        ! ===================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 6
           real(qp),parameter :: hundred = 100.0_qp
           real(qp),parameter :: meighth = -0.125_qp
           real(qp),parameter :: piover2 = 1.57079632679489661923132169163975144210_qp

           ! Local Scalars
           logical(lk) :: colmajor,lquery,restart11,restart12,restart21,restart22,wantu1, &
                     wantu2,wantv1t,wantv2t
           integer(ilp) :: i,imin,imax,iter,iu1cs,iu1sn,iu2cs,iu2sn,iv1tcs,iv1tsn, &
                     iv2tcs,iv2tsn,j,lrworkmin,lrworkopt,maxit,mini
           real(qp) :: b11bulge,b12bulge,b21bulge,b22bulge,dummy,eps,mu,nu,r,sigma11, &
                     sigma21,temp,thetamax,thetamin,thresh,tol,tolmul,unfl,x1,x2,y1,y2
           ! Intrinsic Functions
           intrinsic :: abs,atan2,cos,max,min,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lrwork == -1
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           if (m < 0) then
              info = -6
           else if (p < 0 .or. p > m) then
              info = -7
           else if (q < 0 .or. q > m) then
              info = -8
           else if (q > p .or. q > m - p .or. q > m - q) then
              info = -8
           else if (wantu1 .and. ldu1 < p) then
              info = -12
           else if (wantu2 .and. ldu2 < m - p) then
              info = -14
           else if (wantv1t .and. ldv1t < q) then
              info = -16
           else if (wantv2t .and. ldv2t < m - q) then
              info = -18
           end if
           ! quick return if q = 0
           if (info == 0 .and. q == 0) then
              lrworkmin = 1
              rwork(1) = lrworkmin
              return
           end if
           ! compute workspace
           if (info == 0) then
              iu1cs = 1
              iu1sn = iu1cs + q
              iu2cs = iu1sn + q
              iu2sn = iu2cs + q
              iv1tcs = iu2sn + q
              iv1tsn = iv1tcs + q
              iv2tcs = iv1tsn + q
              iv2tsn = iv2tcs + q
              lrworkopt = iv2tsn + q - 1
              lrworkmin = lrworkopt
              rwork(1) = lrworkopt
              if (lrwork < lrworkmin .and. .not. lquery) then
                 info = -28
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WBBCSD',-info)
              return
           else if (lquery) then
              return
           end if
           ! get machine constants
           eps = la_qlamch('EPSILON')
           unfl = la_qlamch('SAFE MINIMUM')
           tolmul = max(ten,min(hundred,eps**meighth))
           tol = tolmul*eps
           thresh = max(tol,maxitr*q*q*unfl)
           ! test for negligible sines or cosines
           do i = 1,q
              if (theta(i) < thresh) then
                 theta(i) = zero
              else if (theta(i) > piover2 - thresh) then
                 theta(i) = piover2
              end if
           end do
           do i = 1,q - 1
              if (phi(i) < thresh) then
                 phi(i) = zero
              else if (phi(i) > piover2 - thresh) then
                 phi(i) = piover2
              end if
           end do
           ! initial deflation
           imax = q
           do while (imax > 1)
              if (phi(imax - 1) /= zero) then
                 exit
              end if
              imax = imax - 1
           end do
           imin = imax - 1
           if (imin > 1) then
              do while (phi(imin - 1) /= zero)
                 imin = imin - 1
                 if (imin <= 1) exit
              end do
           end if
           ! initialize iteration counter
           maxit = maxitr*q*q
           iter = 0
           ! begin main iteration loop
           do while (imax > 1)
              ! compute the matrix entries
              b11d(imin) = cos(theta(imin))
              b21d(imin) = -sin(theta(imin))
              do i = imin,imax - 1
                 b11e(i) = -sin(theta(i))*sin(phi(i))
                 b11d(i + 1) = cos(theta(i + 1))*cos(phi(i))
                 b12d(i) = sin(theta(i))*cos(phi(i))
                 b12e(i) = cos(theta(i + 1))*sin(phi(i))
                 b21e(i) = -cos(theta(i))*sin(phi(i))
                 b21d(i + 1) = -sin(theta(i + 1))*cos(phi(i))
                 b22d(i) = cos(theta(i))*cos(phi(i))
                 b22e(i) = -sin(theta(i + 1))*sin(phi(i))
              end do
              b12d(imax) = sin(theta(imax))
              b22d(imax) = cos(theta(imax))
              ! abort if not converging; otherwise, increment iter
              if (iter > maxit) then
                 info = 0
                 do i = 1,q
                    if (phi(i) /= zero) info = info + 1
                 end do
                 return
              end if
              iter = iter + imax - imin
              ! compute shifts
              thetamax = theta(imin)
              thetamin = theta(imin)
              do i = imin + 1,imax
                 if (theta(i) > thetamax) thetamax = theta(i)
                 if (theta(i) < thetamin) thetamin = theta(i)
              end do
              if (thetamax > piover2 - thresh) then
                 ! zero on diagonals of b11 and b22; induce deflation with a
                 ! zero shift
                 mu = zero
                 nu = one
              else if (thetamin < thresh) then
                 ! zero on diagonals of b12 and b22; induce deflation with a
                 ! zero shift
                 mu = one
                 nu = zero
              else
                 ! compute shifts for b11 and b21 and use the lesser
                 call la_qlas2(b11d(imax - 1),b11e(imax - 1),b11d(imax),sigma11,dummy)

                 call la_qlas2(b21d(imax - 1),b21e(imax - 1),b21d(imax),sigma21,dummy)

                 if (sigma11 <= sigma21) then
                    mu = sigma11
                    nu = sqrt(one - mu**2)
                    if (mu < thresh) then
                       mu = zero
                       nu = one
                    end if
                 else
                    nu = sigma21
                    mu = sqrt(1.0_qp - nu**2)
                    if (nu < thresh) then
                       mu = one
                       nu = zero
                    end if
                 end if
              end if
              ! rotate to produce bulges in b11 and b21
              if (mu <= nu) then
                 call la_qlartgs(b11d(imin),b11e(imin),mu,rwork(iv1tcs + imin - 1),rwork( &
                           iv1tsn + imin - 1))
              else
                 call la_qlartgs(b21d(imin),b21e(imin),nu,rwork(iv1tcs + imin - 1),rwork( &
                           iv1tsn + imin - 1))
              end if
              temp = rwork(iv1tcs + imin - 1)*b11d(imin) + rwork(iv1tsn + imin - 1)*b11e(imin)
              b11e(imin) = rwork(iv1tcs + imin - 1)*b11e(imin) - rwork(iv1tsn + imin - 1)*b11d(imin)

              b11d(imin) = temp
              b11bulge = rwork(iv1tsn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = rwork(iv1tcs + imin - 1)*b11d(imin + 1)
              temp = rwork(iv1tcs + imin - 1)*b21d(imin) + rwork(iv1tsn + imin - 1)*b21e(imin)
              b21e(imin) = rwork(iv1tcs + imin - 1)*b21e(imin) - rwork(iv1tsn + imin - 1)*b21d(imin)

              b21d(imin) = temp
              b21bulge = rwork(iv1tsn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = rwork(iv1tcs + imin - 1)*b21d(imin + 1)
              ! compute theta(imin)
              theta(imin) = atan2(sqrt(b21d(imin)**2 + b21bulge**2),sqrt(b11d(imin)**2 + &
                        b11bulge**2))
              ! chase the bulges in b11(imin+1,imin) and b21(imin+1,imin)
              if (b11d(imin)**2 + b11bulge**2 > thresh**2) then
                 call la_qlartgp(b11bulge,b11d(imin),rwork(iu1sn + imin - 1),rwork(iu1cs + imin - &
                           1),r)
              else if (mu <= nu) then
                 call la_qlartgs(b11e(imin),b11d(imin + 1),mu,rwork(iu1cs + imin - 1), &
                           rwork(iu1sn + imin - 1))
              else
                 call la_qlartgs(b12d(imin),b12e(imin),nu,rwork(iu1cs + imin - 1),rwork( &
                           iu1sn + imin - 1))
              end if
              if (b21d(imin)**2 + b21bulge**2 > thresh**2) then
                 call la_qlartgp(b21bulge,b21d(imin),rwork(iu2sn + imin - 1),rwork(iu2cs + imin - &
                           1),r)
              else if (nu < mu) then
                 call la_qlartgs(b21e(imin),b21d(imin + 1),nu,rwork(iu2cs + imin - 1), &
                           rwork(iu2sn + imin - 1))
              else
                 call la_qlartgs(b22d(imin),b22e(imin),mu,rwork(iu2cs + imin - 1),rwork(iu2sn + &
                           imin - 1))
              end if
              rwork(iu2cs + imin - 1) = -rwork(iu2cs + imin - 1)
              rwork(iu2sn + imin - 1) = -rwork(iu2sn + imin - 1)
              temp = rwork(iu1cs + imin - 1)*b11e(imin) + rwork(iu1sn + imin - 1)*b11d(imin + 1)
              b11d(imin + 1) = rwork(iu1cs + imin - 1)*b11d(imin + 1) - rwork(iu1sn + imin - 1)*b11e(imin)

              b11e(imin) = temp
              if (imax > imin + 1) then
                 b11bulge = rwork(iu1sn + imin - 1)*b11e(imin + 1)
                 b11e(imin + 1) = rwork(iu1cs + imin - 1)*b11e(imin + 1)
              end if
              temp = rwork(iu1cs + imin - 1)*b12d(imin) + rwork(iu1sn + imin - 1)*b12e(imin)
              b12e(imin) = rwork(iu1cs + imin - 1)*b12e(imin) - rwork(iu1sn + imin - 1)*b12d(imin)
              b12d(imin) = temp
              b12bulge = rwork(iu1sn + imin - 1)*b12d(imin + 1)
              b12d(imin + 1) = rwork(iu1cs + imin - 1)*b12d(imin + 1)
              temp = rwork(iu2cs + imin - 1)*b21e(imin) + rwork(iu2sn + imin - 1)*b21d(imin + 1)
              b21d(imin + 1) = rwork(iu2cs + imin - 1)*b21d(imin + 1) - rwork(iu2sn + imin - 1)*b21e(imin)

              b21e(imin) = temp
              if (imax > imin + 1) then
                 b21bulge = rwork(iu2sn + imin - 1)*b21e(imin + 1)
                 b21e(imin + 1) = rwork(iu2cs + imin - 1)*b21e(imin + 1)
              end if
              temp = rwork(iu2cs + imin - 1)*b22d(imin) + rwork(iu2sn + imin - 1)*b22e(imin)
              b22e(imin) = rwork(iu2cs + imin - 1)*b22e(imin) - rwork(iu2sn + imin - 1)*b22d(imin)
              b22d(imin) = temp
              b22bulge = rwork(iu2sn + imin - 1)*b22d(imin + 1)
              b22d(imin + 1) = rwork(iu2cs + imin - 1)*b22d(imin + 1)
              ! inner loop: chase bulges from b11(imin,imin+2),
              ! b12(imin,imin+1), b21(imin,imin+2), and b22(imin,imin+1) to
              ! bottom-right
              do i = imin + 1,imax - 1
                 ! compute phi(i-1)
                 x1 = sin(theta(i - 1))*b11e(i - 1) + cos(theta(i - 1))*b21e(i - 1)
                 x2 = sin(theta(i - 1))*b11bulge + cos(theta(i - 1))*b21bulge
                 y1 = sin(theta(i - 1))*b12d(i - 1) + cos(theta(i - 1))*b22d(i - 1)
                 y2 = sin(theta(i - 1))*b12bulge + cos(theta(i - 1))*b22bulge
                 phi(i - 1) = atan2(sqrt(x1**2 + x2**2),sqrt(y1**2 + y2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11e(i - 1)**2 + b11bulge**2 <= thresh**2
                 restart21 = b21e(i - 1)**2 + b21bulge**2 <= thresh**2
                 restart12 = b12d(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart22 = b22d(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i-1,i+1), b12(i-1,i),
                 ! b21(i-1,i+1), and b22(i-1,i). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart21) then
                    call la_qlartgp(x2,x1,rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1),r)
                 else if (.not. restart11 .and. restart21) then
                    call la_qlartgp(b11bulge,b11e(i - 1),rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1), &
                               r)
                 else if (restart11 .and. .not. restart21) then
                    call la_qlartgp(b21bulge,b21e(i - 1),rwork(iv1tsn + i - 1),rwork(iv1tcs + i - 1), &
                               r)
                 else if (mu <= nu) then
                    call la_qlartgs(b11d(i),b11e(i),mu,rwork(iv1tcs + i - 1),rwork(iv1tsn + i - 1) &
                               )
                 else
                    call la_qlartgs(b21d(i),b21e(i),nu,rwork(iv1tcs + i - 1),rwork(iv1tsn + i - 1) &
                               )
                 end if
                 rwork(iv1tcs + i - 1) = -rwork(iv1tcs + i - 1)
                 rwork(iv1tsn + i - 1) = -rwork(iv1tsn + i - 1)
                 if (.not. restart12 .and. .not. restart22) then
                    call la_qlartgp(y2,y1,rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - 1 - 1),r)

                 else if (.not. restart12 .and. restart22) then
                    call la_qlartgp(b12bulge,b12d(i - 1),rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - &
                              1 - 1),r)
                 else if (restart12 .and. .not. restart22) then
                    call la_qlartgp(b22bulge,b22d(i - 1),rwork(iv2tsn + i - 1 - 1),rwork(iv2tcs + i - &
                              1 - 1),r)
                 else if (nu < mu) then
                    call la_qlartgs(b12e(i - 1),b12d(i),nu,rwork(iv2tcs + i - 1 - 1),rwork(iv2tsn + &
                              i - 1 - 1))
                 else
                    call la_qlartgs(b22e(i - 1),b22d(i),mu,rwork(iv2tcs + i - 1 - 1),rwork(iv2tsn + &
                              i - 1 - 1))
                 end if
                 temp = rwork(iv1tcs + i - 1)*b11d(i) + rwork(iv1tsn + i - 1)*b11e(i)
                 b11e(i) = rwork(iv1tcs + i - 1)*b11e(i) - rwork(iv1tsn + i - 1)*b11d(i)
                 b11d(i) = temp
                 b11bulge = rwork(iv1tsn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = rwork(iv1tcs + i - 1)*b11d(i + 1)
                 temp = rwork(iv1tcs + i - 1)*b21d(i) + rwork(iv1tsn + i - 1)*b21e(i)
                 b21e(i) = rwork(iv1tcs + i - 1)*b21e(i) - rwork(iv1tsn + i - 1)*b21d(i)
                 b21d(i) = temp
                 b21bulge = rwork(iv1tsn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = rwork(iv1tcs + i - 1)*b21d(i + 1)
                 temp = rwork(iv2tcs + i - 1 - 1)*b12e(i - 1) + rwork(iv2tsn + i - 1 - 1)*b12d(i)
                 b12d(i) = rwork(iv2tcs + i - 1 - 1)*b12d(i) - rwork(iv2tsn + i - 1 - 1)*b12e(i - 1)
                 b12e(i - 1) = temp
                 b12bulge = rwork(iv2tsn + i - 1 - 1)*b12e(i)
                 b12e(i) = rwork(iv2tcs + i - 1 - 1)*b12e(i)
                 temp = rwork(iv2tcs + i - 1 - 1)*b22e(i - 1) + rwork(iv2tsn + i - 1 - 1)*b22d(i)
                 b22d(i) = rwork(iv2tcs + i - 1 - 1)*b22d(i) - rwork(iv2tsn + i - 1 - 1)*b22e(i - 1)
                 b22e(i - 1) = temp
                 b22bulge = rwork(iv2tsn + i - 1 - 1)*b22e(i)
                 b22e(i) = rwork(iv2tcs + i - 1 - 1)*b22e(i)
                 ! compute theta(i)
                 x1 = cos(phi(i - 1))*b11d(i) + sin(phi(i - 1))*b12e(i - 1)
                 x2 = cos(phi(i - 1))*b11bulge + sin(phi(i - 1))*b12bulge
                 y1 = cos(phi(i - 1))*b21d(i) + sin(phi(i - 1))*b22e(i - 1)
                 y2 = cos(phi(i - 1))*b21bulge + sin(phi(i - 1))*b22bulge
                 theta(i) = atan2(sqrt(y1**2 + y2**2),sqrt(x1**2 + x2**2))
                 ! determine if there are bulges to chase or if a new direct
                 ! summand has been reached
                 restart11 = b11d(i)**2 + b11bulge**2 <= thresh**2
                 restart12 = b12e(i - 1)**2 + b12bulge**2 <= thresh**2
                 restart21 = b21d(i)**2 + b21bulge**2 <= thresh**2
                 restart22 = b22e(i - 1)**2 + b22bulge**2 <= thresh**2
                 ! if possible, chase bulges from b11(i+1,i), b12(i+1,i-1),
                 ! b21(i+1,i), and b22(i+1,i-1). if necessary, restart bulge-
                 ! chasing by applying the original shift again.
                 if (.not. restart11 .and. .not. restart12) then
                    call la_qlartgp(x2,x1,rwork(iu1sn + i - 1),rwork(iu1cs + i - 1),r)
                 else if (.not. restart11 .and. restart12) then
                    call la_qlartgp(b11bulge,b11d(i),rwork(iu1sn + i - 1),rwork(iu1cs + i - 1),r)

                 else if (restart11 .and. .not. restart12) then
                    call la_qlartgp(b12bulge,b12e(i - 1),rwork(iu1sn + i - 1),rwork(iu1cs + i - 1), &
                              r)
                 else if (mu <= nu) then
                    call la_qlartgs(b11e(i),b11d(i + 1),mu,rwork(iu1cs + i - 1),rwork(iu1sn + i - 1) &
                               )
                 else
                    call la_qlartgs(b12d(i),b12e(i),nu,rwork(iu1cs + i - 1),rwork(iu1sn + i - 1))

                 end if
                 if (.not. restart21 .and. .not. restart22) then
                    call la_qlartgp(y2,y1,rwork(iu2sn + i - 1),rwork(iu2cs + i - 1),r)
                 else if (.not. restart21 .and. restart22) then
                    call la_qlartgp(b21bulge,b21d(i),rwork(iu2sn + i - 1),rwork(iu2cs + i - 1),r)

                 else if (restart21 .and. .not. restart22) then
                    call la_qlartgp(b22bulge,b22e(i - 1),rwork(iu2sn + i - 1),rwork(iu2cs + i - 1), &
                              r)
                 else if (nu < mu) then
                    call la_qlartgs(b21e(i),b21e(i + 1),nu,rwork(iu2cs + i - 1),rwork(iu2sn + i - 1) &
                               )
                 else
                    call la_qlartgs(b22d(i),b22e(i),mu,rwork(iu2cs + i - 1),rwork(iu2sn + i - 1))

                 end if
                 rwork(iu2cs + i - 1) = -rwork(iu2cs + i - 1)
                 rwork(iu2sn + i - 1) = -rwork(iu2sn + i - 1)
                 temp = rwork(iu1cs + i - 1)*b11e(i) + rwork(iu1sn + i - 1)*b11d(i + 1)
                 b11d(i + 1) = rwork(iu1cs + i - 1)*b11d(i + 1) - rwork(iu1sn + i - 1)*b11e(i)
                 b11e(i) = temp
                 if (i < imax - 1) then
                    b11bulge = rwork(iu1sn + i - 1)*b11e(i + 1)
                    b11e(i + 1) = rwork(iu1cs + i - 1)*b11e(i + 1)
                 end if
                 temp = rwork(iu2cs + i - 1)*b21e(i) + rwork(iu2sn + i - 1)*b21d(i + 1)
                 b21d(i + 1) = rwork(iu2cs + i - 1)*b21d(i + 1) - rwork(iu2sn + i - 1)*b21e(i)
                 b21e(i) = temp
                 if (i < imax - 1) then
                    b21bulge = rwork(iu2sn + i - 1)*b21e(i + 1)
                    b21e(i + 1) = rwork(iu2cs + i - 1)*b21e(i + 1)
                 end if
                 temp = rwork(iu1cs + i - 1)*b12d(i) + rwork(iu1sn + i - 1)*b12e(i)
                 b12e(i) = rwork(iu1cs + i - 1)*b12e(i) - rwork(iu1sn + i - 1)*b12d(i)
                 b12d(i) = temp
                 b12bulge = rwork(iu1sn + i - 1)*b12d(i + 1)
                 b12d(i + 1) = rwork(iu1cs + i - 1)*b12d(i + 1)
                 temp = rwork(iu2cs + i - 1)*b22d(i) + rwork(iu2sn + i - 1)*b22e(i)
                 b22e(i) = rwork(iu2cs + i - 1)*b22e(i) - rwork(iu2sn + i - 1)*b22d(i)
                 b22d(i) = temp
                 b22bulge = rwork(iu2sn + i - 1)*b22d(i + 1)
                 b22d(i + 1) = rwork(iu2cs + i - 1)*b22d(i + 1)
              end do
              ! compute phi(imax-1)
              x1 = sin(theta(imax - 1))*b11e(imax - 1) + cos(theta(imax - 1))*b21e(imax - 1)
              y1 = sin(theta(imax - 1))*b12d(imax - 1) + cos(theta(imax - 1))*b22d(imax - 1)
              y2 = sin(theta(imax - 1))*b12bulge + cos(theta(imax - 1))*b22bulge
              phi(imax - 1) = atan2(abs(x1),sqrt(y1**2 + y2**2))
              ! chase bulges from b12(imax-1,imax) and b22(imax-1,imax)
              restart12 = b12d(imax - 1)**2 + b12bulge**2 <= thresh**2
              restart22 = b22d(imax - 1)**2 + b22bulge**2 <= thresh**2
              if (.not. restart12 .and. .not. restart22) then
                 call la_qlartgp(y2,y1,rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + imax - 1 - 1),r)

              else if (.not. restart12 .and. restart22) then
                 call la_qlartgp(b12bulge,b12d(imax - 1),rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + &
                           imax - 1 - 1),r)
              else if (restart12 .and. .not. restart22) then
                 call la_qlartgp(b22bulge,b22d(imax - 1),rwork(iv2tsn + imax - 1 - 1),rwork(iv2tcs + &
                           imax - 1 - 1),r)
              else if (nu < mu) then
                 call la_qlartgs(b12e(imax - 1),b12d(imax),nu,rwork(iv2tcs + imax - 1 - 1),rwork( &
                           iv2tsn + imax - 1 - 1))
              else
                 call la_qlartgs(b22e(imax - 1),b22d(imax),mu,rwork(iv2tcs + imax - 1 - 1),rwork( &
                           iv2tsn + imax - 1 - 1))
              end if
              temp = rwork(iv2tcs + imax - 1 - 1)*b12e(imax - 1) + rwork(iv2tsn + imax - 1 - 1)*b12d(imax)

              b12d(imax) = rwork(iv2tcs + imax - 1 - 1)*b12d(imax) - rwork(iv2tsn + imax - 1 - 1)*b12e(imax - 1)

              b12e(imax - 1) = temp
              temp = rwork(iv2tcs + imax - 1 - 1)*b22e(imax - 1) + rwork(iv2tsn + imax - 1 - 1)*b22d(imax)

              b22d(imax) = rwork(iv2tcs + imax - 1 - 1)*b22d(imax) - rwork(iv2tsn + imax - 1 - 1)*b22e(imax - 1)

              b22e(imax - 1) = temp
              ! update singular vectors
              if (wantu1) then
                 if (colmajor) then
                    call la_wlasr('R','V','F',p,imax - imin + 1,rwork(iu1cs + imin - 1),rwork( &
                              iu1sn + imin - 1),u1(1,imin),ldu1)
                 else
                    call la_wlasr('L','V','F',imax - imin + 1,p,rwork(iu1cs + imin - 1),rwork( &
                              iu1sn + imin - 1),u1(imin,1),ldu1)
                 end if
              end if
              if (wantu2) then
                 if (colmajor) then
                    call la_wlasr('R','V','F',m - p,imax - imin + 1,rwork(iu2cs + imin - 1),rwork( &
                              iu2sn + imin - 1),u2(1,imin),ldu2)
                 else
                    call la_wlasr('L','V','F',imax - imin + 1,m - p,rwork(iu2cs + imin - 1),rwork( &
                              iu2sn + imin - 1),u2(imin,1),ldu2)
                 end if
              end if
              if (wantv1t) then
                 if (colmajor) then
                    call la_wlasr('L','V','F',imax - imin + 1,q,rwork(iv1tcs + imin - 1),rwork( &
                              iv1tsn + imin - 1),v1t(imin,1),ldv1t)
                 else
                    call la_wlasr('R','V','F',q,imax - imin + 1,rwork(iv1tcs + imin - 1),rwork( &
                              iv1tsn + imin - 1),v1t(1,imin),ldv1t)
                 end if
              end if
              if (wantv2t) then
                 if (colmajor) then
                    call la_wlasr('L','V','F',imax - imin + 1,m - q,rwork(iv2tcs + imin - 1), &
                              rwork(iv2tsn + imin - 1),v2t(imin,1),ldv2t)
                 else
                    call la_wlasr('R','V','F',m - q,imax - imin + 1,rwork(iv2tcs + imin - 1), &
                              rwork(iv2tsn + imin - 1),v2t(1,imin),ldv2t)
                 end if
              end if
              ! fix signs on b11(imax-1,imax) and b21(imax-1,imax)
              if (b11e(imax - 1) + b21e(imax - 1) > 0) then
                 b11d(imax) = -b11d(imax)
                 b21d(imax) = -b21d(imax)
                 if (wantv1t) then
                    if (colmajor) then
                       call la_wscal(q,cnegone,v1t(imax,1),ldv1t)
                    else
                       call la_wscal(q,cnegone,v1t(1,imax),1)
                    end if
                 end if
              end if
              ! compute theta(imax)
              x1 = cos(phi(imax - 1))*b11d(imax) + sin(phi(imax - 1))*b12e(imax - 1)
              y1 = cos(phi(imax - 1))*b21d(imax) + sin(phi(imax - 1))*b22e(imax - 1)
              theta(imax) = atan2(abs(y1),abs(x1))
              ! fix signs on b11(imax,imax), b12(imax,imax-1), b21(imax,imax),
              ! and b22(imax,imax-1)
              if (b11d(imax) + b12e(imax - 1) < 0) then
                 b12d(imax) = -b12d(imax)
                 if (wantu1) then
                    if (colmajor) then
                       call la_wscal(p,cnegone,u1(1,imax),1)
                    else
                       call la_wscal(p,cnegone,u1(imax,1),ldu1)
                    end if
                 end if
              end if
              if (b21d(imax) + b22e(imax - 1) > 0) then
                 b22d(imax) = -b22d(imax)
                 if (wantu2) then
                    if (colmajor) then
                       call la_wscal(m - p,cnegone,u2(1,imax),1)
                    else
                       call la_wscal(m - p,cnegone,u2(imax,1),ldu2)
                    end if
                 end if
              end if
              ! fix signs on b12(imax,imax) and b22(imax,imax)
              if (b12d(imax) + b22d(imax) < 0) then
                 if (wantv2t) then
                    if (colmajor) then
                       call la_wscal(m - q,cnegone,v2t(imax,1),ldv2t)
                    else
                       call la_wscal(m - q,cnegone,v2t(1,imax),1)
                    end if
                 end if
              end if
              ! test for negligible sines or cosines
              do i = imin,imax
                 if (theta(i) < thresh) then
                    theta(i) = zero
                 else if (theta(i) > piover2 - thresh) then
                    theta(i) = piover2
                 end if
              end do
              do i = imin,imax - 1
                 if (phi(i) < thresh) then
                    phi(i) = zero
                 else if (phi(i) > piover2 - thresh) then
                    phi(i) = piover2
                 end if
              end do
              ! deflate
              if (imax > 1) then
                 do while (phi(imax - 1) == zero)
                    imax = imax - 1
                    if (imax <= 1) exit
                 end do
              end if
              if (imin > imax - 1) imin = imax - 1
              if (imin > 1) then
                 do while (phi(imin - 1) /= zero)
                     imin = imin - 1
                     if (imin <= 1) exit
                 end do
              end if
              ! repeat main iteration loop
           end do
           ! postprocessing: order theta from least to greatest
           do i = 1,q
              mini = i
              thetamin = theta(i)
              do j = i + 1,q
                 if (theta(j) < thetamin) then
                    mini = j
                    thetamin = theta(j)
                 end if
              end do
              if (mini /= i) then
                 theta(mini) = theta(i)
                 theta(i) = thetamin
                 if (colmajor) then
                    if (wantu1) call la_wswap(p,u1(1,i),1,u1(1,mini),1)
                    if (wantu2) call la_wswap(m - p,u2(1,i),1,u2(1,mini),1)
                    if (wantv1t) call la_wswap(q,v1t(i,1),ldv1t,v1t(mini,1),ldv1t)

                    if (wantv2t) call la_wswap(m - q,v2t(i,1),ldv2t,v2t(mini,1),ldv2t)

                 else
                    if (wantu1) call la_wswap(p,u1(i,1),ldu1,u1(mini,1),ldu1)
                    if (wantu2) call la_wswap(m - p,u2(i,1),ldu2,u2(mini,1),ldu2)
                    if (wantv1t) call la_wswap(q,v1t(1,i),1,v1t(1,mini),1)
                    if (wantv2t) call la_wswap(m - q,v2t(1,i),1,v2t(1,mini),1)
                 end if
              end if
           end do
           return
     end subroutine la_wbbcsd

     !> CUNBDB5: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then some other vector from the orthogonal complement
     !> is returned. This vector is chosen in an arbitrary but deterministic
     !> way.

     pure subroutine la_cunbdb5(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           complex(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,j
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('CUNBDB5',-info)
              return
           end if
           ! project x onto the orthogonal complement of q
           call la_cunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work,lwork, &
                     childinfo)
           ! if the projection is nonzero, then return
           if (la_scnrm2(m1,x1,incx1) /= czero .or. la_scnrm2(m2,x2,incx2) /= czero) &
                     then
              return
           end if
           ! project each standard basis vector e_1,...,e_m1 in turn, stopping
           ! when a nonzero projection is found
           do i = 1,m1
              do j = 1,m1
                 x1(j) = czero
              end do
              x1(i) = cone
              do j = 1,m2
                 x2(j) = czero
              end do
              call la_cunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_scnrm2(m1,x1,incx1) /= czero .or. la_scnrm2(m2,x2,incx2) /= czero) &
                        then
                 return
              end if
           end do
           ! project each standard basis vector e_(m1+1),...,e_(m1+m2) in turn,
           ! stopping when a nonzero projection is found
           do i = 1,m2
              do j = 1,m1
                 x1(j) = czero
              end do
              do j = 1,m2
                 x2(j) = czero
              end do
              x2(i) = cone
              call la_cunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_scnrm2(m1,x1,incx1) /= czero .or. la_scnrm2(m2,x2,incx2) /= czero) &
                        then
                 return
              end if
           end do
           return
     end subroutine la_cunbdb5
     !> ZUNBDB5: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then some other vector from the orthogonal complement
     !> is returned. This vector is chosen in an arbitrary but deterministic
     !> way.

     pure subroutine la_zunbdb5(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           complex(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,j
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('ZUNBDB5',-info)
              return
           end if
           ! project x onto the orthogonal complement of q
           call la_zunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work,lwork, &
                     childinfo)
           ! if the projection is nonzero, then return
           if (la_dznrm2(m1,x1,incx1) /= czero .or. la_dznrm2(m2,x2,incx2) /= czero) &
                     then
              return
           end if
           ! project each standard basis vector e_1,...,e_m1 in turn, stopping
           ! when a nonzero projection is found
           do i = 1,m1
              do j = 1,m1
                 x1(j) = czero
              end do
              x1(i) = cone
              do j = 1,m2
                 x2(j) = czero
              end do
              call la_zunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_dznrm2(m1,x1,incx1) /= czero .or. la_dznrm2(m2,x2,incx2) /= czero) &
                        then
                 return
              end if
           end do
           ! project each standard basis vector e_(m1+1),...,e_(m1+m2) in turn,
           ! stopping when a nonzero projection is found
           do i = 1,m2
              do j = 1,m1
                 x1(j) = czero
              end do
              do j = 1,m2
                 x2(j) = czero
              end do
              x2(i) = cone
              call la_zunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_dznrm2(m1,x1,incx1) /= czero .or. la_dznrm2(m2,x2,incx2) /= czero) &
                        then
                 return
              end if
           end do
           return
     end subroutine la_zunbdb5
     !> WUNBDB5: orthogonalizes the column vector
     !> X = [ X1 ]
     !> [ X2 ]
     !> with respect to the columns of
     !> Q = [ Q1 ] .
     !> [ Q2 ]
     !> The columns of Q must be orthonormal.
     !> If the projection is zero according to Kahan's "twice is enough"
     !> criterion, then some other vector from the orthogonal complement
     !> is returned. This vector is chosen in an arbitrary but deterministic
     !> way.

     pure subroutine la_wunbdb5(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
               lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx1,incx2,ldq1,ldq2,lwork,m1,m2,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(in) :: q1(ldq1,*),q2(ldq2,*)
           complex(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: x1(*),x2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,j
           ! Intrinsic Function
           intrinsic :: max
           ! Executable Statements
           ! test input arguments
           info = 0
           if (m1 < 0) then
              info = -1
           else if (m2 < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (incx1 < 1) then
              info = -5
           else if (incx2 < 1) then
              info = -7
           else if (ldq1 < max(1,m1)) then
              info = -9
           else if (ldq2 < max(1,m2)) then
              info = -11
           else if (lwork < n) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('WUNBDB5',-info)
              return
           end if
           ! project x onto the orthogonal complement of q
           call la_wunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work,lwork, &
                     childinfo)
           ! if the projection is nonzero, then return
           if (la_qwnrm2(m1,x1,incx1) /= czero .or. la_qwnrm2(m2,x2,incx2) /= czero) &
                     then
              return
           end if
           ! project each standard basis vector e_1,...,e_m1 in turn, stopping
           ! when a nonzero projection is found
           do i = 1,m1
              do j = 1,m1
                 x1(j) = czero
              end do
              x1(i) = cone
              do j = 1,m2
                 x2(j) = czero
              end do
              call la_wunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_qwnrm2(m1,x1,incx1) /= czero .or. la_qwnrm2(m2,x2,incx2) /= czero) &
                        then
                 return
              end if
           end do
           ! project each standard basis vector e_(m1+1),...,e_(m1+m2) in turn,
           ! stopping when a nonzero projection is found
           do i = 1,m2
              do j = 1,m1
                 x1(j) = czero
              end do
              do j = 1,m2
                 x2(j) = czero
              end do
              x2(i) = cone
              call la_wunbdb6(m1,m2,n,x1,incx1,x2,incx2,q1,ldq1,q2,ldq2,work, &
                        lwork,childinfo)
              if (la_qwnrm2(m1,x1,incx1) /= czero .or. la_qwnrm2(m2,x2,incx2) /= czero) &
                        then
                 return
              end if
           end do
           return
     end subroutine la_wunbdb5

     !> CUNCSD: computes the CS decomposition of an M-by-M partitioned
     !> unitary matrix X:
     !> [  I  0  0 |  0  0  0 ]
     !> [  0  C  0 |  0 -S  0 ]
     !> [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
     !> X = [-----------] = [---------] [---------------------] [---------]   .
     !> [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
     !> [  0  S  0 |  0  C  0 ]
     !> [  0  0  I |  0  0  0 ]
     !> X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P,
     !> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     !> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     !> which R = MIN(P,M-P,Q,M-Q).

     recursive subroutine la_cuncsd(jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,p,q,x11, &
     ldx11,x12,ldx12,x21,ldx21,x22,ldx22,theta,u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t, &
               work,lwork,rwork,lrwork,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,ldx11,ldx12,ldx21,ldx22, &
                     lrwork,lwork,m,p,q
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: theta(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*),work(*)

           complex(sp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ===================================================================

           ! Local Scalars
           character :: transt,signst
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,itauq2,j,lbbcsdwork, &
           lbbcsdworkmin,lbbcsdworkopt,lorbdbwork,lorbdbworkmin,lorbdbworkopt,lorglqwork, &
           lorglqworkmin,lorglqworkopt,lorgqrwork,lorgqrworkmin,lorgqrworkopt,lworkmin, &
                     lworkopt,p1,q1
           logical(lk) :: colmajor,defaultsigns,lquery,wantu1,wantu2,wantv1t,wantv2t
           integer(ilp) :: lrworkmin,lrworkopt
           logical(lk) :: lrquery
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           defaultsigns = .not. la_lsame(signs,'O')
           lquery = lwork == -1
           lrquery = lrwork == -1
           if (m < 0) then
              info = -7
           else if (p < 0 .or. p > m) then
              info = -8
           else if (q < 0 .or. q > m) then
              info = -9
           else if (colmajor .and. ldx11 < max(1,p)) then
             info = -11
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
             info = -11
           else if (colmajor .and. ldx12 < max(1,p)) then
             info = -13
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
             info = -13
           else if (colmajor .and. ldx21 < max(1,m - p)) then
             info = -15
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
             info = -15
           else if (colmajor .and. ldx22 < max(1,m - p)) then
             info = -17
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
             info = -17
           else if (wantu1 .and. ldu1 < p) then
              info = -20
           else if (wantu2 .and. ldu2 < m - p) then
              info = -22
           else if (wantv1t .and. ldv1t < q) then
              info = -24
           else if (wantv2t .and. ldv2t < m - q) then
              info = -26
           end if
           ! work with transpose if convenient
           if (info == 0 .and. min(p,m - p) < min(q,m - q)) then
              if (colmajor) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_cuncsd(jobv1t,jobv2t,jobu1,jobu2,transt,signst,m,q,p,x11, &
              ldx11,x21,ldx21,x12,ldx12,x22,ldx22,theta,v1t,ldv1t,v2t,ldv2t,u1,ldu1, &
                        u2,ldu2,work,lwork,rwork,lrwork,iwork,info)
              return
           end if
           ! work with permutation [ 0 i; i 0 ] * x * [ 0 i; i 0 ] if
           ! convenient
           if (info == 0 .and. m - q < q) then
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_cuncsd(jobu2,jobu1,jobv2t,jobv1t,trans,signst,m,m - p,m - q,x22, &
              ldx22,x21,ldx21,x12,ldx12,x11,ldx11,theta,u2,ldu2,u1,ldu1,v2t,ldv2t, &
                        v1t,ldv1t,work,lwork,rwork,lrwork,iwork,info)
              return
           end if
           ! compute workspace
           if (info == 0) then
              ! real workspace
              iphi = 2
              ib11d = iphi + max(1,q - 1)
              ib11e = ib11d + max(1,q)
              ib12d = ib11e + max(1,q - 1)
              ib12e = ib12d + max(1,q)
              ib21d = ib12e + max(1,q - 1)
              ib21e = ib21d + max(1,q)
              ib22d = ib21e + max(1,q - 1)
              ib22e = ib22d + max(1,q)
              ibbcsd = ib22e + max(1,q - 1)
              call la_cbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,theta,u1, &
              ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,theta,theta,theta,theta,theta,theta, &
                        theta,theta,rwork,-1,childinfo)
              lbbcsdworkopt = int(rwork(1),KIND=ilp)
              lbbcsdworkmin = lbbcsdworkopt
              lrworkopt = ibbcsd + lbbcsdworkopt - 1
              lrworkmin = ibbcsd + lbbcsdworkmin - 1
              rwork(1) = lrworkopt
              ! complex workspace
              itaup1 = 2
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              itauq2 = itauq1 + max(1,q)
              iorgqr = itauq2 + max(1,m - q)
              call la_cungqr(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorgqrworkopt = int(work(1),KIND=ilp)
              lorgqrworkmin = max(1,m - q)
              iorglq = itauq2 + max(1,m - q)
              call la_cunglq(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorglqworkopt = int(work(1),KIND=ilp)
              lorglqworkmin = max(1,m - q)
              iorbdb = itauq2 + max(1,m - q)
              call la_cunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
                        ldx22,theta,theta,u1,u2,v1t,v2t,work,-1,childinfo)
              lorbdbworkopt = int(work(1),KIND=ilp)
              lorbdbworkmin = lorbdbworkopt
              lworkopt = max(iorgqr + lorgqrworkopt,iorglq + lorglqworkopt,iorbdb + &
                        lorbdbworkopt) - 1
              lworkmin = max(iorgqr + lorgqrworkmin,iorglq + lorglqworkmin,iorbdb + &
                        lorbdbworkmin) - 1
              work(1) = max(lworkopt,lworkmin)
              if (lwork < lworkmin .and. .not. (lquery .or. lrquery)) then
                 info = -22
              else if (lrwork < lrworkmin .and. .not. (lquery .or. lrquery)) then
                 info = -24
              else
                 lorgqrwork = lwork - iorgqr + 1
                 lorglqwork = lwork - iorglq + 1
                 lorbdbwork = lwork - iorbdb + 1
                 lbbcsdwork = lrwork - ibbcsd + 1
              end if
           end if
           ! abort if any illegal arguments
           if (info /= 0) then
              call la_xerbla('CUNCSD',-info)
              return
           else if (lquery .or. lrquery) then
              return
           end if
           ! transform to bidiagonal block form
           call la_cunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
           ldx22,theta,rwork(iphi),work(itaup1),work(itaup2),work(itauq1),work(itauq2),work( &
                     iorbdb),lorbdbwork,childinfo)
           ! accumulate householder reflectors
           if (colmajor) then
              if (wantu1 .and. p > 0) then
                 call la_clacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_cungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqrwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_clacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_cungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqrwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_clacpy('U',q - 1,q - 1,x11(1,2),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_cunglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglqwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_clacpy('U',p,m - q,x12,ldx12,v2t,ldv2t)
                 if (m - p > q) then
                    call la_clacpy('U',m - p - q,m - p - q,x22(q + 1,p + 1),ldx22,v2t(p + 1,p + 1), &
                              ldv2t)
                 end if
                 if (m > q) then
                    call la_cunglq(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorglq), &
                              lorglqwork,info)
                 end if
              end if
           else
              if (wantu1 .and. p > 0) then
                 call la_clacpy('U',q,p,x11,ldx11,u1,ldu1)
                 call la_cunglq(p,p,q,u1,ldu1,work(itaup1),work(iorglq),lorglqwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_clacpy('U',q,m - p,x21,ldx21,u2,ldu2)
                 call la_cunglq(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorglq),lorglqwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_clacpy('L',q - 1,q - 1,x11(2,1),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_cungqr(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorgqr), &
                           lorgqrwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 p1 = min(p + 1,m)
                 q1 = min(q + 1,m)
                 call la_clacpy('L',m - q,p,x12,ldx12,v2t,ldv2t)
                 if (m > p + q) then
                    call la_clacpy('L',m - p - q,m - p - q,x22(p1,q1),ldx22,v2t(p + 1,p + 1),ldv2t)

                 end if
                 call la_cungqr(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorgqr), &
                           lorgqrwork,info)
              end if
           end if
           ! compute the csd of the matrix in bidiagonal-block form
           call la_cbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,rwork(iphi), &
           u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
           rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                     lbbcsdwork,info)
           ! permute rows and columns to place identity submatrices in top-
           ! left corner of (1,1)-block and/or bottom-right corner of (1,2)-
           ! block and/or bottom-right corner of (2,1)-block and/or top-left
           ! corner of (2,2)-block
           if (q > 0 .and. wantu2) then
              do i = 1,q
                 iwork(i) = m - p - q + i
              end do
              do i = q + 1,m - p
                 iwork(i) = i - q
              end do
              if (colmajor) then
                 call la_clapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              else
                 call la_clapmr(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           end if
           if (m > 0 .and. wantv2t) then
              do i = 1,p
                 iwork(i) = m - p - q + i
              end do
              do i = p + 1,m - q
                 iwork(i) = i - p
              end do
              if (.not. colmajor) then
                 call la_clapmt(.false.,m - q,m - q,v2t,ldv2t,iwork)
              else
                 call la_clapmr(.false.,m - q,m - q,v2t,ldv2t,iwork)
              end if
           end if
           return
           ! end la_cuncsd
     end subroutine la_cuncsd
     !> ZUNCSD: computes the CS decomposition of an M-by-M partitioned
     !> unitary matrix X:
     !> [  I  0  0 |  0  0  0 ]
     !> [  0  C  0 |  0 -S  0 ]
     !> [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
     !> X = [-----------] = [---------] [---------------------] [---------]   .
     !> [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
     !> [  0  S  0 |  0  C  0 ]
     !> [  0  0  I |  0  0  0 ]
     !> X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P,
     !> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     !> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     !> which R = MIN(P,M-P,Q,M-Q).

     recursive subroutine la_zuncsd(jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,p,q,x11, &
     ldx11,x12,ldx12,x21,ldx21,x22,ldx22,theta,u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t, &
               work,lwork,rwork,lrwork,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,ldx11,ldx12,ldx21,ldx22, &
                     lrwork,lwork,m,p,q
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: theta(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*),work(*)

           complex(dp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ===================================================================

           ! Local Scalars
           character :: transt,signst
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,itauq2,j,lbbcsdwork, &
           lbbcsdworkmin,lbbcsdworkopt,lorbdbwork,lorbdbworkmin,lorbdbworkopt,lorglqwork, &
           lorglqworkmin,lorglqworkopt,lorgqrwork,lorgqrworkmin,lorgqrworkopt,lworkmin, &
                     lworkopt,p1,q1
           logical(lk) :: colmajor,defaultsigns,lquery,wantu1,wantu2,wantv1t,wantv2t
           integer(ilp) :: lrworkmin,lrworkopt
           logical(lk) :: lrquery
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           defaultsigns = .not. la_lsame(signs,'O')
           lquery = lwork == -1
           lrquery = lrwork == -1
           if (m < 0) then
              info = -7
           else if (p < 0 .or. p > m) then
              info = -8
           else if (q < 0 .or. q > m) then
              info = -9
           else if (colmajor .and. ldx11 < max(1,p)) then
             info = -11
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
             info = -11
           else if (colmajor .and. ldx12 < max(1,p)) then
             info = -13
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
             info = -13
           else if (colmajor .and. ldx21 < max(1,m - p)) then
             info = -15
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
             info = -15
           else if (colmajor .and. ldx22 < max(1,m - p)) then
             info = -17
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
             info = -17
           else if (wantu1 .and. ldu1 < p) then
              info = -20
           else if (wantu2 .and. ldu2 < m - p) then
              info = -22
           else if (wantv1t .and. ldv1t < q) then
              info = -24
           else if (wantv2t .and. ldv2t < m - q) then
              info = -26
           end if
           ! work with transpose if convenient
           if (info == 0 .and. min(p,m - p) < min(q,m - q)) then
              if (colmajor) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_zuncsd(jobv1t,jobv2t,jobu1,jobu2,transt,signst,m,q,p,x11, &
              ldx11,x21,ldx21,x12,ldx12,x22,ldx22,theta,v1t,ldv1t,v2t,ldv2t,u1,ldu1, &
                        u2,ldu2,work,lwork,rwork,lrwork,iwork,info)
              return
           end if
           ! work with permutation [ 0 i; i 0 ] * x * [ 0 i; i 0 ] if
           ! convenient
           if (info == 0 .and. m - q < q) then
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_zuncsd(jobu2,jobu1,jobv2t,jobv1t,trans,signst,m,m - p,m - q,x22, &
              ldx22,x21,ldx21,x12,ldx12,x11,ldx11,theta,u2,ldu2,u1,ldu1,v2t,ldv2t, &
                        v1t,ldv1t,work,lwork,rwork,lrwork,iwork,info)
              return
           end if
           ! compute workspace
           if (info == 0) then
              ! real workspace
              iphi = 2
              ib11d = iphi + max(1,q - 1)
              ib11e = ib11d + max(1,q)
              ib12d = ib11e + max(1,q - 1)
              ib12e = ib12d + max(1,q)
              ib21d = ib12e + max(1,q - 1)
              ib21e = ib21d + max(1,q)
              ib22d = ib21e + max(1,q - 1)
              ib22e = ib22d + max(1,q)
              ibbcsd = ib22e + max(1,q - 1)
              call la_zbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,theta,u1, &
              ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,theta,theta,theta,theta,theta,theta, &
                        theta,theta,rwork,-1,childinfo)
              lbbcsdworkopt = int(rwork(1),KIND=ilp)
              lbbcsdworkmin = lbbcsdworkopt
              lrworkopt = ibbcsd + lbbcsdworkopt - 1
              lrworkmin = ibbcsd + lbbcsdworkmin - 1
              rwork(1) = lrworkopt
              ! complex workspace
              itaup1 = 2
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              itauq2 = itauq1 + max(1,q)
              iorgqr = itauq2 + max(1,m - q)
              call la_zungqr(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorgqrworkopt = int(work(1),KIND=ilp)
              lorgqrworkmin = max(1,m - q)
              iorglq = itauq2 + max(1,m - q)
              call la_zunglq(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorglqworkopt = int(work(1),KIND=ilp)
              lorglqworkmin = max(1,m - q)
              iorbdb = itauq2 + max(1,m - q)
              call la_zunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
                        ldx22,theta,theta,u1,u2,v1t,v2t,work,-1,childinfo)
              lorbdbworkopt = int(work(1),KIND=ilp)
              lorbdbworkmin = lorbdbworkopt
              lworkopt = max(iorgqr + lorgqrworkopt,iorglq + lorglqworkopt,iorbdb + &
                        lorbdbworkopt) - 1
              lworkmin = max(iorgqr + lorgqrworkmin,iorglq + lorglqworkmin,iorbdb + &
                        lorbdbworkmin) - 1
              work(1) = max(lworkopt,lworkmin)
              if (lwork < lworkmin .and. .not. (lquery .or. lrquery)) then
                 info = -22
              else if (lrwork < lrworkmin .and. .not. (lquery .or. lrquery)) then
                 info = -24
              else
                 lorgqrwork = lwork - iorgqr + 1
                 lorglqwork = lwork - iorglq + 1
                 lorbdbwork = lwork - iorbdb + 1
                 lbbcsdwork = lrwork - ibbcsd + 1
              end if
           end if
           ! abort if any illegal arguments
           if (info /= 0) then
              call la_xerbla('ZUNCSD',-info)
              return
           else if (lquery .or. lrquery) then
              return
           end if
           ! transform to bidiagonal block form
           call la_zunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
           ldx22,theta,rwork(iphi),work(itaup1),work(itaup2),work(itauq1),work(itauq2),work( &
                     iorbdb),lorbdbwork,childinfo)
           ! accumulate householder reflectors
           if (colmajor) then
              if (wantu1 .and. p > 0) then
                 call la_zlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_zungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqrwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_zlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_zungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqrwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_zlacpy('U',q - 1,q - 1,x11(1,2),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_zunglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglqwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_zlacpy('U',p,m - q,x12,ldx12,v2t,ldv2t)
                 if (m - p > q) then
                    call la_zlacpy('U',m - p - q,m - p - q,x22(q + 1,p + 1),ldx22,v2t(p + 1,p + 1), &
                              ldv2t)
                 end if
                 if (m > q) then
                    call la_zunglq(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorglq), &
                              lorglqwork,info)
                 end if
              end if
           else
              if (wantu1 .and. p > 0) then
                 call la_zlacpy('U',q,p,x11,ldx11,u1,ldu1)
                 call la_zunglq(p,p,q,u1,ldu1,work(itaup1),work(iorglq),lorglqwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_zlacpy('U',q,m - p,x21,ldx21,u2,ldu2)
                 call la_zunglq(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorglq),lorglqwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_zlacpy('L',q - 1,q - 1,x11(2,1),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_zungqr(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorgqr), &
                           lorgqrwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 p1 = min(p + 1,m)
                 q1 = min(q + 1,m)
                 call la_zlacpy('L',m - q,p,x12,ldx12,v2t,ldv2t)
                 if (m > p + q) then
                    call la_zlacpy('L',m - p - q,m - p - q,x22(p1,q1),ldx22,v2t(p + 1,p + 1),ldv2t)

                 end if
                 call la_zungqr(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorgqr), &
                           lorgqrwork,info)
              end if
           end if
           ! compute the csd of the matrix in bidiagonal-block form
           call la_zbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,rwork(iphi), &
           u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
           rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                     lbbcsdwork,info)
           ! permute rows and columns to place identity submatrices in top-
           ! left corner of (1,1)-block and/or bottom-right corner of (1,2)-
           ! block and/or bottom-right corner of (2,1)-block and/or top-left
           ! corner of (2,2)-block
           if (q > 0 .and. wantu2) then
              do i = 1,q
                 iwork(i) = m - p - q + i
              end do
              do i = q + 1,m - p
                 iwork(i) = i - q
              end do
              if (colmajor) then
                 call la_zlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              else
                 call la_zlapmr(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           end if
           if (m > 0 .and. wantv2t) then
              do i = 1,p
                 iwork(i) = m - p - q + i
              end do
              do i = p + 1,m - q
                 iwork(i) = i - p
              end do
              if (.not. colmajor) then
                 call la_zlapmt(.false.,m - q,m - q,v2t,ldv2t,iwork)
              else
                 call la_zlapmr(.false.,m - q,m - q,v2t,ldv2t,iwork)
              end if
           end if
           return
           ! end la_zuncsd
     end subroutine la_zuncsd
     !> WUNCSD: computes the CS decomposition of an M-by-M partitioned
     !> unitary matrix X:
     !> [  I  0  0 |  0  0  0 ]
     !> [  0  C  0 |  0 -S  0 ]
     !> [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
     !> X = [-----------] = [---------] [---------------------] [---------]   .
     !> [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
     !> [  0  S  0 |  0  C  0 ]
     !> [  0  0  I |  0  0  0 ]
     !> X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P,
     !> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
     !> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
     !> which R = MIN(P,M-P,Q,M-Q).

     recursive subroutine la_wuncsd(jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,p,q,x11, &
     ldx11,x12,ldx12,x21,ldx21,x22,ldx22,theta,u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t, &
               work,lwork,rwork,lrwork,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t,jobv2t,signs,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,ldv2t,ldx11,ldx12,ldx21,ldx22, &
                     lrwork,lwork,m,p,q
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: theta(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),v2t(ldv2t,*),work(*)

           complex(qp),intent(inout) :: x11(ldx11,*),x12(ldx12,*),x21(ldx21,*),x22(ldx22,*)

        ! ===================================================================

           ! Local Scalars
           character :: transt,signst
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,itauq2,j,lbbcsdwork, &
           lbbcsdworkmin,lbbcsdworkopt,lorbdbwork,lorbdbworkmin,lorbdbworkopt,lorglqwork, &
           lorglqworkmin,lorglqworkopt,lorgqrwork,lorgqrworkmin,lorgqrworkopt,lworkmin, &
                     lworkopt,p1,q1
           logical(lk) :: colmajor,defaultsigns,lquery,wantu1,wantu2,wantv1t,wantv2t
           integer(ilp) :: lrworkmin,lrworkopt
           logical(lk) :: lrquery
           ! Intrinsic Functions
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           wantv2t = la_lsame(jobv2t,'Y')
           colmajor = .not. la_lsame(trans,'T')
           defaultsigns = .not. la_lsame(signs,'O')
           lquery = lwork == -1
           lrquery = lrwork == -1
           if (m < 0) then
              info = -7
           else if (p < 0 .or. p > m) then
              info = -8
           else if (q < 0 .or. q > m) then
              info = -9
           else if (colmajor .and. ldx11 < max(1,p)) then
             info = -11
           else if (.not. colmajor .and. ldx11 < max(1,q)) then
             info = -11
           else if (colmajor .and. ldx12 < max(1,p)) then
             info = -13
           else if (.not. colmajor .and. ldx12 < max(1,m - q)) then
             info = -13
           else if (colmajor .and. ldx21 < max(1,m - p)) then
             info = -15
           else if (.not. colmajor .and. ldx21 < max(1,q)) then
             info = -15
           else if (colmajor .and. ldx22 < max(1,m - p)) then
             info = -17
           else if (.not. colmajor .and. ldx22 < max(1,m - q)) then
             info = -17
           else if (wantu1 .and. ldu1 < p) then
              info = -20
           else if (wantu2 .and. ldu2 < m - p) then
              info = -22
           else if (wantv1t .and. ldv1t < q) then
              info = -24
           else if (wantv2t .and. ldv2t < m - q) then
              info = -26
           end if
           ! work with transpose if convenient
           if (info == 0 .and. min(p,m - p) < min(q,m - q)) then
              if (colmajor) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_wuncsd(jobv1t,jobv2t,jobu1,jobu2,transt,signst,m,q,p,x11, &
              ldx11,x21,ldx21,x12,ldx12,x22,ldx22,theta,v1t,ldv1t,v2t,ldv2t,u1,ldu1, &
                        u2,ldu2,work,lwork,rwork,lrwork,iwork,info)
              return
           end if
           ! work with permutation [ 0 i; i 0 ] * x * [ 0 i; i 0 ] if
           ! convenient
           if (info == 0 .and. m - q < q) then
              if (defaultsigns) then
                 signst = 'O'
              else
                 signst = 'D'
              end if
              call la_wuncsd(jobu2,jobu1,jobv2t,jobv1t,trans,signst,m,m - p,m - q,x22, &
              ldx22,x21,ldx21,x12,ldx12,x11,ldx11,theta,u2,ldu2,u1,ldu1,v2t,ldv2t, &
                        v1t,ldv1t,work,lwork,rwork,lrwork,iwork,info)
              return
           end if
           ! compute workspace
           if (info == 0) then
              ! real workspace
              iphi = 2
              ib11d = iphi + max(1,q - 1)
              ib11e = ib11d + max(1,q)
              ib12d = ib11e + max(1,q - 1)
              ib12e = ib12d + max(1,q)
              ib21d = ib12e + max(1,q - 1)
              ib21e = ib21d + max(1,q)
              ib22d = ib21e + max(1,q - 1)
              ib22e = ib22d + max(1,q)
              ibbcsd = ib22e + max(1,q - 1)
              call la_wbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,theta,u1, &
              ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,theta,theta,theta,theta,theta,theta, &
                        theta,theta,rwork,-1,childinfo)
              lbbcsdworkopt = int(rwork(1),KIND=ilp)
              lbbcsdworkmin = lbbcsdworkopt
              lrworkopt = ibbcsd + lbbcsdworkopt - 1
              lrworkmin = ibbcsd + lbbcsdworkmin - 1
              rwork(1) = lrworkopt
              ! complex workspace
              itaup1 = 2
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              itauq2 = itauq1 + max(1,q)
              iorgqr = itauq2 + max(1,m - q)
              call la_wungqr(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorgqrworkopt = int(work(1),KIND=ilp)
              lorgqrworkmin = max(1,m - q)
              iorglq = itauq2 + max(1,m - q)
              call la_wunglq(m - q,m - q,m - q,u1,max(1,m - q),u1,work,-1,childinfo)
              lorglqworkopt = int(work(1),KIND=ilp)
              lorglqworkmin = max(1,m - q)
              iorbdb = itauq2 + max(1,m - q)
              call la_wunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
                        ldx22,theta,theta,u1,u2,v1t,v2t,work,-1,childinfo)
              lorbdbworkopt = int(work(1),KIND=ilp)
              lorbdbworkmin = lorbdbworkopt
              lworkopt = max(iorgqr + lorgqrworkopt,iorglq + lorglqworkopt,iorbdb + &
                        lorbdbworkopt) - 1
              lworkmin = max(iorgqr + lorgqrworkmin,iorglq + lorglqworkmin,iorbdb + &
                        lorbdbworkmin) - 1
              work(1) = max(lworkopt,lworkmin)
              if (lwork < lworkmin .and. .not. (lquery .or. lrquery)) then
                 info = -22
              else if (lrwork < lrworkmin .and. .not. (lquery .or. lrquery)) then
                 info = -24
              else
                 lorgqrwork = lwork - iorgqr + 1
                 lorglqwork = lwork - iorglq + 1
                 lorbdbwork = lwork - iorbdb + 1
                 lbbcsdwork = lrwork - ibbcsd + 1
              end if
           end if
           ! abort if any illegal arguments
           if (info /= 0) then
              call la_xerbla('WUNCSD',-info)
              return
           else if (lquery .or. lrquery) then
              return
           end if
           ! transform to bidiagonal block form
           call la_wunbdb(trans,signs,m,p,q,x11,ldx11,x12,ldx12,x21,ldx21,x22, &
           ldx22,theta,rwork(iphi),work(itaup1),work(itaup2),work(itauq1),work(itauq2),work( &
                     iorbdb),lorbdbwork,childinfo)
           ! accumulate householder reflectors
           if (colmajor) then
              if (wantu1 .and. p > 0) then
                 call la_wlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_wungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqrwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_wlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_wungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqrwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_wlacpy('U',q - 1,q - 1,x11(1,2),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_wunglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglqwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 call la_wlacpy('U',p,m - q,x12,ldx12,v2t,ldv2t)
                 if (m - p > q) then
                    call la_wlacpy('U',m - p - q,m - p - q,x22(q + 1,p + 1),ldx22,v2t(p + 1,p + 1), &
                              ldv2t)
                 end if
                 if (m > q) then
                    call la_wunglq(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorglq), &
                              lorglqwork,info)
                 end if
              end if
           else
              if (wantu1 .and. p > 0) then
                 call la_wlacpy('U',q,p,x11,ldx11,u1,ldu1)
                 call la_wunglq(p,p,q,u1,ldu1,work(itaup1),work(iorglq),lorglqwork, &
                           info)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_wlacpy('U',q,m - p,x21,ldx21,u2,ldu2)
                 call la_wunglq(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorglq),lorglqwork, &
                            info)
              end if
              if (wantv1t .and. q > 0) then
                 call la_wlacpy('L',q - 1,q - 1,x11(2,1),ldx11,v1t(2,2),ldv1t)
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_wungqr(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorgqr), &
                           lorgqrwork,info)
              end if
              if (wantv2t .and. m - q > 0) then
                 p1 = min(p + 1,m)
                 q1 = min(q + 1,m)
                 call la_wlacpy('L',m - q,p,x12,ldx12,v2t,ldv2t)
                 if (m > p + q) then
                    call la_wlacpy('L',m - p - q,m - p - q,x22(p1,q1),ldx22,v2t(p + 1,p + 1),ldv2t)

                 end if
                 call la_wungqr(m - q,m - q,m - q,v2t,ldv2t,work(itauq2),work(iorgqr), &
                           lorgqrwork,info)
              end if
           end if
           ! compute the csd of the matrix in bidiagonal-block form
           call la_wbbcsd(jobu1,jobu2,jobv1t,jobv2t,trans,m,p,q,theta,rwork(iphi), &
           u1,ldu1,u2,ldu2,v1t,ldv1t,v2t,ldv2t,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
           rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                     lbbcsdwork,info)
           ! permute rows and columns to place identity submatrices in top-
           ! left corner of (1,1)-block and/or bottom-right corner of (1,2)-
           ! block and/or bottom-right corner of (2,1)-block and/or top-left
           ! corner of (2,2)-block
           if (q > 0 .and. wantu2) then
              do i = 1,q
                 iwork(i) = m - p - q + i
              end do
              do i = q + 1,m - p
                 iwork(i) = i - q
              end do
              if (colmajor) then
                 call la_wlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              else
                 call la_wlapmr(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           end if
           if (m > 0 .and. wantv2t) then
              do i = 1,p
                 iwork(i) = m - p - q + i
              end do
              do i = p + 1,m - q
                 iwork(i) = i - p
              end do
              if (.not. colmajor) then
                 call la_wlapmt(.false.,m - q,m - q,v2t,ldv2t,iwork)
              else
                 call la_wlapmr(.false.,m - q,m - q,v2t,ldv2t,iwork)
              end if
           end if
           return
           ! end la_wuncsd
     end subroutine la_wuncsd

     !> CUNBDB1: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
     !> M-P, or M-Q. Routines CUNBDB2, CUNBDB3, and CUNBDB4 handle cases in
     !> which Q is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_cunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           complex(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < q .or. m - p < q) then
              info = -2
           else if (q < 0 .or. m - q < q) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 2
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CUNBDB1',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., q of x11 and x21
           do i = 1,q
              call la_clarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              call la_clarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              theta(i) = atan2(real(x21(i,i),KIND=sp),real(x11(i,i),KIND=sp))
              c = cos(theta(i))
              s = sin(theta(i))
              x11(i,i) = cone
              x21(i,i) = cone
              call la_clarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
              call la_clarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
              if (i < q) then
                 call la_csrot(q - i,x11(i,i + 1),ldx11,x21(i,i + 1),ldx21,c,s)
                 call la_clacgv(q - i,x21(i,i + 1),ldx21)
                 call la_clarfgp(q - i,x21(i,i + 1),x21(i,i + 2),ldx21,tauq1(i))
                 s = real(x21(i,i + 1),KIND=sp)
                 x21(i,i + 1) = cone
                 call la_clarf('R',p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
                 call la_clarf('R',m - p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
                 call la_clacgv(q - i,x21(i,i + 1),ldx21)
                 c = sqrt(la_scnrm2(p - i,x11(i + 1,i + 1),1)**2 + la_scnrm2(m - p - i,x21(i + &
                           1,i + 1),1)**2)
                 phi(i) = atan2(s,c)
                 call la_cunbdb5(p - i,m - p - i,q - i - 1,x11(i + 1,i + 1),1,x21(i + 1,i + 1),1,x11(i + 1, &
                           i + 2),ldx11,x21(i + 1,i + 2),ldx21,work(iorbdb5),lorbdb5,childinfo)
              end if
           end do
           return
     end subroutine la_cunbdb1
     !> ZUNBDB1: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
     !> M-P, or M-Q. Routines ZUNBDB2, ZUNBDB3, and ZUNBDB4 handle cases in
     !> which Q is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_zunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           complex(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < q .or. m - p < q) then
              info = -2
           else if (q < 0 .or. m - q < q) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 2
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZUNBDB1',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., q of x11 and x21
           do i = 1,q
              call la_zlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              theta(i) = atan2(real(x21(i,i),KIND=dp),real(x11(i,i),KIND=dp))
              c = cos(theta(i))
              s = sin(theta(i))
              x11(i,i) = cone
              x21(i,i) = cone
              call la_zlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
              call la_zlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
              if (i < q) then
                 call la_zdrot(q - i,x11(i,i + 1),ldx11,x21(i,i + 1),ldx21,c,s)
                 call la_zlacgv(q - i,x21(i,i + 1),ldx21)
                 call la_zlarfgp(q - i,x21(i,i + 1),x21(i,i + 2),ldx21,tauq1(i))
                 s = real(x21(i,i + 1),KIND=dp)
                 x21(i,i + 1) = cone
                 call la_zlarf('R',p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
                 call la_zlarf('R',m - p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
                 call la_zlacgv(q - i,x21(i,i + 1),ldx21)
                 c = sqrt(la_dznrm2(p - i,x11(i + 1,i + 1),1)**2 + la_dznrm2(m - p - i,x21(i + &
                           1,i + 1),1)**2)
                 phi(i) = atan2(s,c)
                 call la_zunbdb5(p - i,m - p - i,q - i - 1,x11(i + 1,i + 1),1,x21(i + 1,i + 1),1,x11(i + 1, &
                           i + 2),ldx11,x21(i + 1,i + 2),ldx21,work(iorbdb5),lorbdb5,childinfo)
              end if
           end do
           return
     end subroutine la_zunbdb1
     !> WUNBDB1: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
     !> M-P, or M-Q. Routines WUNBDB2, WUNBDB3, and WUNBDB4 handle cases in
     !> which Q is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_wunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           complex(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < q .or. m - p < q) then
              info = -2
           else if (q < 0 .or. m - q < q) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 2
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WUNBDB1',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., q of x11 and x21
           do i = 1,q
              call la_wlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              theta(i) = atan2(real(x21(i,i),KIND=qp),real(x11(i,i),KIND=qp))
              c = cos(theta(i))
              s = sin(theta(i))
              x11(i,i) = cone
              x21(i,i) = cone
              call la_wlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
              call la_wlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
              if (i < q) then
                 call la_wqrot(q - i,x11(i,i + 1),ldx11,x21(i,i + 1),ldx21,c,s)
                 call la_wlacgv(q - i,x21(i,i + 1),ldx21)
                 call la_wlarfgp(q - i,x21(i,i + 1),x21(i,i + 2),ldx21,tauq1(i))
                 s = real(x21(i,i + 1),KIND=qp)
                 x21(i,i + 1) = cone
                 call la_wlarf('R',p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
                 call la_wlarf('R',m - p - i,q - i,x21(i,i + 1),ldx21,tauq1(i),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
                 call la_wlacgv(q - i,x21(i,i + 1),ldx21)
                 c = sqrt(la_qwnrm2(p - i,x11(i + 1,i + 1),1)**2 + la_qwnrm2(m - p - i,x21(i + &
                           1,i + 1),1)**2)
                 phi(i) = atan2(s,c)
                 call la_wunbdb5(p - i,m - p - i,q - i - 1,x11(i + 1,i + 1),1,x21(i + 1,i + 1),1,x11(i + 1, &
                           i + 2),ldx11,x21(i + 1,i + 2),ldx21,work(iorbdb5),lorbdb5,childinfo)
              end if
           end do
           return
     end subroutine la_wunbdb1

     !> CUNBDB2: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
     !> Q, or M-Q. Routines CUNBDB1, CUNBDB3, and CUNBDB4 handle cases in
     !> which P is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_cunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           complex(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < 0 .or. p > m - p) then
              info = -2
           else if (q < 0 .or. q < p .or. m - q < p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CUNBDB2',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., p of x11 and x21
           do i = 1,p
              if (i > 1) then
                 call la_csrot(q - i + 1,x11(i,i),ldx11,x21(i - 1,i),ldx21,c,s)
              end if
              call la_clacgv(q - i + 1,x11(i,i),ldx11)
              call la_clarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              c = real(x11(i,i),KIND=sp)
              x11(i,i) = cone
              call la_clarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_clarf('R',m - p - i + 1,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(i,i),ldx21, &
                        work(ilarf))
              call la_clacgv(q - i + 1,x11(i,i),ldx11)
              s = sqrt(la_scnrm2(p - i,x11(i + 1,i),1)**2 + la_scnrm2(m - p - i + 1,x21(i,i), &
                        1)**2)
              theta(i) = atan2(s,c)
              call la_cunbdb5(p - i,m - p - i + 1,q - i,x11(i + 1,i),1,x21(i,i),1,x11(i + 1,i + 1), &
                        ldx11,x21(i,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_cscal(p - i,cnegone,x11(i + 1,i),1)
              call la_clarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              if (i < p) then
                 call la_clarfgp(p - i,x11(i + 1,i),x11(i + 2,i),1,taup1(i))
                 phi(i) = atan2(real(x11(i + 1,i),KIND=sp),real(x21(i,i),KIND=sp))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x11(i + 1,i) = cone
                 call la_clarf('L',p - i,q - i,x11(i + 1,i),1,conjg(taup1(i)),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
              end if
              x21(i,i) = cone
              call la_clarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
           end do
           ! reduce the bottom-right portion of x21 to the identity matrix
           do i = p + 1,q
              call la_clarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              x21(i,i) = cone
              call la_clarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
           end do
           return
     end subroutine la_cunbdb2
     !> ZUNBDB2: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
     !> Q, or M-Q. Routines ZUNBDB1, ZUNBDB3, and ZUNBDB4 handle cases in
     !> which P is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_zunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           complex(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < 0 .or. p > m - p) then
              info = -2
           else if (q < 0 .or. q < p .or. m - q < p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZUNBDB2',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., p of x11 and x21
           do i = 1,p
              if (i > 1) then
                 call la_zdrot(q - i + 1,x11(i,i),ldx11,x21(i - 1,i),ldx21,c,s)
              end if
              call la_zlacgv(q - i + 1,x11(i,i),ldx11)
              call la_zlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              c = real(x11(i,i),KIND=dp)
              x11(i,i) = cone
              call la_zlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_zlarf('R',m - p - i + 1,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(i,i),ldx21, &
                        work(ilarf))
              call la_zlacgv(q - i + 1,x11(i,i),ldx11)
              s = sqrt(la_dznrm2(p - i,x11(i + 1,i),1)**2 + la_dznrm2(m - p - i + 1,x21(i,i), &
                        1)**2)
              theta(i) = atan2(s,c)
              call la_zunbdb5(p - i,m - p - i + 1,q - i,x11(i + 1,i),1,x21(i,i),1,x11(i + 1,i + 1), &
                        ldx11,x21(i,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_zscal(p - i,cnegone,x11(i + 1,i),1)
              call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              if (i < p) then
                 call la_zlarfgp(p - i,x11(i + 1,i),x11(i + 2,i),1,taup1(i))
                 phi(i) = atan2(real(x11(i + 1,i),KIND=dp),real(x21(i,i),KIND=dp))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x11(i + 1,i) = cone
                 call la_zlarf('L',p - i,q - i,x11(i + 1,i),1,conjg(taup1(i)),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
              end if
              x21(i,i) = cone
              call la_zlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
           end do
           ! reduce the bottom-right portion of x21 to the identity matrix
           do i = p + 1,q
              call la_zlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              x21(i,i) = cone
              call la_zlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
           end do
           return
     end subroutine la_zunbdb2
     !> WUNBDB2: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
     !> Q, or M-Q. Routines WUNBDB1, WUNBDB3, and WUNBDB4 handle cases in
     !> which P is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
     !> angles THETA, PHI.

     subroutine la_wunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           complex(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < 0 .or. p > m - p) then
              info = -2
           else if (q < 0 .or. q < p .or. m - q < p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p - 1,m - p,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WUNBDB2',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., p of x11 and x21
           do i = 1,p
              if (i > 1) then
                 call la_wqrot(q - i + 1,x11(i,i),ldx11,x21(i - 1,i),ldx21,c,s)
              end if
              call la_wlacgv(q - i + 1,x11(i,i),ldx11)
              call la_wlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              c = real(x11(i,i),KIND=qp)
              x11(i,i) = cone
              call la_wlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_wlarf('R',m - p - i + 1,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(i,i),ldx21, &
                        work(ilarf))
              call la_wlacgv(q - i + 1,x11(i,i),ldx11)
              s = sqrt(la_qwnrm2(p - i,x11(i + 1,i),1)**2 + la_qwnrm2(m - p - i + 1,x21(i,i), &
                        1)**2)
              theta(i) = atan2(s,c)
              call la_wunbdb5(p - i,m - p - i + 1,q - i,x11(i + 1,i),1,x21(i,i),1,x11(i + 1,i + 1), &
                        ldx11,x21(i,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_wscal(p - i,cnegone,x11(i + 1,i),1)
              call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              if (i < p) then
                 call la_wlarfgp(p - i,x11(i + 1,i),x11(i + 2,i),1,taup1(i))
                 phi(i) = atan2(real(x11(i + 1,i),KIND=qp),real(x21(i,i),KIND=qp))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x11(i + 1,i) = cone
                 call la_wlarf('L',p - i,q - i,x11(i + 1,i),1,conjg(taup1(i)),x11(i + 1,i + 1), &
                           ldx11,work(ilarf))
              end if
              x21(i,i) = cone
              call la_wlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
           end do
           ! reduce the bottom-right portion of x21 to the identity matrix
           do i = p + 1,q
              call la_wlarfgp(m - p - i + 1,x21(i,i),x21(i + 1,i),1,taup2(i))
              x21(i,i) = cone
              call la_wlarf('L',m - p - i + 1,q - i,x21(i,i),1,conjg(taup2(i)),x21(i,i + 1), &
                        ldx21,work(ilarf))
           end do
           return
     end subroutine la_wunbdb2

     !> CUNBDB3: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
     !> Q, or M-Q. Routines CUNBDB1, CUNBDB2, and CUNBDB4 handle cases in
     !> which M-P is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_cunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           complex(sp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (2*p < m .or. p > m) then
              info = -2
           else if (q < m - p .or. m - q < m - p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CUNBDB3',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., m-p of x11 and x21
           do i = 1,m - p
              if (i > 1) then
                 call la_csrot(q - i + 1,x11(i - 1,i),ldx11,x21(i,i),ldx11,c,s)
              end if
              call la_clacgv(q - i + 1,x21(i,i),ldx21)
              call la_clarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              s = real(x21(i,i),KIND=sp)
              x21(i,i) = cone
              call la_clarf('R',p - i + 1,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i,i),ldx11, &
                        work(ilarf))
              call la_clarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              call la_clacgv(q - i + 1,x21(i,i),ldx21)
              c = sqrt(la_scnrm2(p - i + 1,x11(i,i),1)**2 + la_scnrm2(m - p - i,x21(i + 1,i), &
                        1)**2)
              theta(i) = atan2(s,c)
              call la_cunbdb5(p - i + 1,m - p - i,q - i,x11(i,i),1,x21(i + 1,i),1,x11(i,i + 1), &
                        ldx11,x21(i + 1,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_clarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              if (i < m - p) then
                 call la_clarfgp(m - p - i,x21(i + 1,i),x21(i + 2,i),1,taup2(i))
                 phi(i) = atan2(real(x21(i + 1,i),KIND=sp),real(x11(i,i),KIND=sp))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x21(i + 1,i) = cone
                 call la_clarf('L',m - p - i,q - i,x21(i + 1,i),1,conjg(taup2(i)),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
              end if
              x11(i,i) = cone
              call la_clarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
           end do
           ! reduce the bottom-right portion of x11 to the identity matrix
           do i = m - p + 1,q
              call la_clarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              x11(i,i) = cone
              call la_clarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
           end do
           return
     end subroutine la_cunbdb3
     !> ZUNBDB3: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
     !> Q, or M-Q. Routines ZUNBDB1, ZUNBDB2, and ZUNBDB4 handle cases in
     !> which M-P is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_zunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           complex(dp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (2*p < m .or. p > m) then
              info = -2
           else if (q < m - p .or. m - q < m - p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZUNBDB3',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., m-p of x11 and x21
           do i = 1,m - p
              if (i > 1) then
                 call la_zdrot(q - i + 1,x11(i - 1,i),ldx11,x21(i,i),ldx11,c,s)
              end if
              call la_zlacgv(q - i + 1,x21(i,i),ldx21)
              call la_zlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              s = real(x21(i,i),KIND=dp)
              x21(i,i) = cone
              call la_zlarf('R',p - i + 1,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i,i),ldx11, &
                        work(ilarf))
              call la_zlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              call la_zlacgv(q - i + 1,x21(i,i),ldx21)
              c = sqrt(la_dznrm2(p - i + 1,x11(i,i),1)**2 + la_dznrm2(m - p - i,x21(i + 1,i), &
                        1)**2)
              theta(i) = atan2(s,c)
              call la_zunbdb5(p - i + 1,m - p - i,q - i,x11(i,i),1,x21(i + 1,i),1,x11(i,i + 1), &
                        ldx11,x21(i + 1,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_zlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              if (i < m - p) then
                 call la_zlarfgp(m - p - i,x21(i + 1,i),x21(i + 2,i),1,taup2(i))
                 phi(i) = atan2(real(x21(i + 1,i),KIND=dp),real(x11(i,i),KIND=dp))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x21(i + 1,i) = cone
                 call la_zlarf('L',m - p - i,q - i,x21(i + 1,i),1,conjg(taup2(i)),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
              end if
              x11(i,i) = cone
              call la_zlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
           end do
           ! reduce the bottom-right portion of x11 to the identity matrix
           do i = m - p + 1,q
              call la_zlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              x11(i,i) = cone
              call la_zlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
           end do
           return
     end subroutine la_zunbdb3
     !> WUNBDB3: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,
     !> Q, or M-Q. Routines WUNBDB1, WUNBDB2, and WUNBDB4 handle cases in
     !> which M-P is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_wunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           complex(qp),intent(out) :: taup1(*),taup2(*),tauq1(*),work(*)
           complex(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (2*p < m .or. p > m) then
              info = -2
           else if (q < m - p .or. m - q < m - p) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(p,m - p - 1,q - 1)
              iorbdb5 = 2
              lorbdb5 = q - 1
              lworkopt = max(ilarf + llarf - 1,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WUNBDB3',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce rows 1, ..., m-p of x11 and x21
           do i = 1,m - p
              if (i > 1) then
                 call la_wqrot(q - i + 1,x11(i - 1,i),ldx11,x21(i,i),ldx11,c,s)
              end if
              call la_wlacgv(q - i + 1,x21(i,i),ldx21)
              call la_wlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              s = real(x21(i,i),KIND=qp)
              x21(i,i) = cone
              call la_wlarf('R',p - i + 1,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i,i),ldx11, &
                        work(ilarf))
              call la_wlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              call la_wlacgv(q - i + 1,x21(i,i),ldx21)
              c = sqrt(la_qwnrm2(p - i + 1,x11(i,i),1)**2 + la_qwnrm2(m - p - i,x21(i + 1,i), &
                        1)**2)
              theta(i) = atan2(s,c)
              call la_wunbdb5(p - i + 1,m - p - i,q - i,x11(i,i),1,x21(i + 1,i),1,x11(i,i + 1), &
                        ldx11,x21(i + 1,i + 1),ldx21,work(iorbdb5),lorbdb5,childinfo)
              call la_wlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              if (i < m - p) then
                 call la_wlarfgp(m - p - i,x21(i + 1,i),x21(i + 2,i),1,taup2(i))
                 phi(i) = atan2(real(x21(i + 1,i),KIND=qp),real(x11(i,i),KIND=qp))
                 c = cos(phi(i))
                 s = sin(phi(i))
                 x21(i + 1,i) = cone
                 call la_wlarf('L',m - p - i,q - i,x21(i + 1,i),1,conjg(taup2(i)),x21(i + 1,i + 1), &
                           ldx21,work(ilarf))
              end if
              x11(i,i) = cone
              call la_wlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
           end do
           ! reduce the bottom-right portion of x11 to the identity matrix
           do i = m - p + 1,q
              call la_wlarfgp(p - i + 1,x11(i,i),x11(i + 1,i),1,taup1(i))
              x11(i,i) = cone
              call la_wlarf('L',p - i + 1,q - i,x11(i,i),1,conjg(taup1(i)),x11(i,i + 1),ldx11, &
                        work(ilarf))
           end do
           return
     end subroutine la_wunbdb3

     !> CUNBDB4: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
     !> M-P, or Q. Routines CUNBDB1, CUNBDB2, and CUNBDB3 handle cases in
     !> which M-Q is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_cunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               phantom,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(sp),intent(out) :: phi(*),theta(*)
           complex(sp),intent(out) :: phantom(*),taup1(*),taup2(*),tauq1(*),work(*)
           complex(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(sp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,j,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < m - q .or. m - p < m - q) then
              info = -2
           else if (q < m - q .or. q > m) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(q - 1,p - 1,m - p - 1)
              iorbdb5 = 2
              lorbdb5 = q
              lworkopt = ilarf + llarf - 1
              lworkopt = max(lworkopt,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CUNBDB4',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., m-q of x11 and x21
           do i = 1,m - q
              if (i == 1) then
                 do j = 1,m
                    phantom(j) = czero
                 end do
                 call la_cunbdb5(p,m - p,q,phantom(1),1,phantom(p + 1),1,x11,ldx11,x21, &
                           ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_cscal(p,cnegone,phantom(1),1)
                 call la_clarfgp(p,phantom(1),phantom(2),1,taup1(1))
                 call la_clarfgp(m - p,phantom(p + 1),phantom(p + 2),1,taup2(1))
                 theta(i) = atan2(real(phantom(1),KIND=sp),real(phantom(p + 1),KIND=sp))

                 c = cos(theta(i))
                 s = sin(theta(i))
                 phantom(1) = cone
                 phantom(p + 1) = cone
                 call la_clarf('L',p,q,phantom(1),1,conjg(taup1(1)),x11,ldx11,work( &
                           ilarf))
                 call la_clarf('L',m - p,q,phantom(p + 1),1,conjg(taup2(1)),x21,ldx21, &
                           work(ilarf))
              else
                 call la_cunbdb5(p - i + 1,m - p - i + 1,q - i + 1,x11(i,i - 1),1,x21(i,i - 1),1,x11(i,i) &
                           ,ldx11,x21(i,i),ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_cscal(p - i + 1,cnegone,x11(i,i - 1),1)
                 call la_clarfgp(p - i + 1,x11(i,i - 1),x11(i + 1,i - 1),1,taup1(i))
                 call la_clarfgp(m - p - i + 1,x21(i,i - 1),x21(i + 1,i - 1),1,taup2(i))
                 theta(i) = atan2(real(x11(i,i - 1),KIND=sp),real(x21(i,i - 1),KIND=sp))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 x11(i,i - 1) = cone
                 x21(i,i - 1) = cone
                 call la_clarf('L',p - i + 1,q - i + 1,x11(i,i - 1),1,conjg(taup1(i)),x11(i,i), &
                           ldx11,work(ilarf))
                 call la_clarf('L',m - p - i + 1,q - i + 1,x21(i,i - 1),1,conjg(taup2(i)),x21(i,i), &
                           ldx21,work(ilarf))
              end if
              call la_csrot(q - i + 1,x11(i,i),ldx11,x21(i,i),ldx21,s,-c)
              call la_clacgv(q - i + 1,x21(i,i),ldx21)
              call la_clarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              c = real(x21(i,i),KIND=sp)
              x21(i,i) = cone
              call la_clarf('R',p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_clarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              call la_clacgv(q - i + 1,x21(i,i),ldx21)
              if (i < m - q) then
                 s = sqrt(la_scnrm2(p - i,x11(i + 1,i),1)**2 + la_scnrm2(m - p - i,x21(i + 1, &
                           i),1)**2)
                 phi(i) = atan2(s,c)
              end if
           end do
           ! reduce the bottom-right portion of x11 to [ i 0 ]
           do i = m - q + 1,p
              call la_clacgv(q - i + 1,x11(i,i),ldx11)
              call la_clarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              x11(i,i) = cone
              call la_clarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_clarf('R',q - p,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(m - q + 1,i),ldx21, &
                        work(ilarf))
              call la_clacgv(q - i + 1,x11(i,i),ldx11)
           end do
           ! reduce the bottom-right portion of x21 to [ 0 i ]
           do i = p + 1,q
              call la_clacgv(q - i + 1,x21(m - q + i - p,i),ldx21)
              call la_clarfgp(q - i + 1,x21(m - q + i - p,i),x21(m - q + i - p,i + 1),ldx21,tauq1(i))

              x21(m - q + i - p,i) = cone
              call la_clarf('R',q - i,q - i + 1,x21(m - q + i - p,i),ldx21,tauq1(i),x21(m - q + i - p + 1,i) &
                        ,ldx21,work(ilarf))
              call la_clacgv(q - i + 1,x21(m - q + i - p,i),ldx21)
           end do
           return
     end subroutine la_cunbdb4
     !> ZUNBDB4: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
     !> M-P, or Q. Routines ZUNBDB1, ZUNBDB2, and ZUNBDB3 handle cases in
     !> which M-Q is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_zunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               phantom,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(dp),intent(out) :: phi(*),theta(*)
           complex(dp),intent(out) :: phantom(*),taup1(*),taup2(*),tauq1(*),work(*)
           complex(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(dp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,j,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < m - q .or. m - p < m - q) then
              info = -2
           else if (q < m - q .or. q > m) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(q - 1,p - 1,m - p - 1)
              iorbdb5 = 2
              lorbdb5 = q
              lworkopt = ilarf + llarf - 1
              lworkopt = max(lworkopt,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZUNBDB4',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., m-q of x11 and x21
           do i = 1,m - q
              if (i == 1) then
                 do j = 1,m
                    phantom(j) = czero
                 end do
                 call la_zunbdb5(p,m - p,q,phantom(1),1,phantom(p + 1),1,x11,ldx11,x21, &
                           ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_zscal(p,cnegone,phantom(1),1)
                 call la_zlarfgp(p,phantom(1),phantom(2),1,taup1(1))
                 call la_zlarfgp(m - p,phantom(p + 1),phantom(p + 2),1,taup2(1))
                 theta(i) = atan2(real(phantom(1),KIND=dp),real(phantom(p + 1),KIND=dp))

                 c = cos(theta(i))
                 s = sin(theta(i))
                 phantom(1) = cone
                 phantom(p + 1) = cone
                 call la_zlarf('L',p,q,phantom(1),1,conjg(taup1(1)),x11,ldx11,work( &
                           ilarf))
                 call la_zlarf('L',m - p,q,phantom(p + 1),1,conjg(taup2(1)),x21,ldx21, &
                           work(ilarf))
              else
                 call la_zunbdb5(p - i + 1,m - p - i + 1,q - i + 1,x11(i,i - 1),1,x21(i,i - 1),1,x11(i,i) &
                           ,ldx11,x21(i,i),ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_zscal(p - i + 1,cnegone,x11(i,i - 1),1)
                 call la_zlarfgp(p - i + 1,x11(i,i - 1),x11(i + 1,i - 1),1,taup1(i))
                 call la_zlarfgp(m - p - i + 1,x21(i,i - 1),x21(i + 1,i - 1),1,taup2(i))
                 theta(i) = atan2(real(x11(i,i - 1),KIND=dp),real(x21(i,i - 1),KIND=dp))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 x11(i,i - 1) = cone
                 x21(i,i - 1) = cone
                 call la_zlarf('L',p - i + 1,q - i + 1,x11(i,i - 1),1,conjg(taup1(i)),x11(i,i), &
                           ldx11,work(ilarf))
                 call la_zlarf('L',m - p - i + 1,q - i + 1,x21(i,i - 1),1,conjg(taup2(i)),x21(i,i), &
                           ldx21,work(ilarf))
              end if
              call la_zdrot(q - i + 1,x11(i,i),ldx11,x21(i,i),ldx21,s,-c)
              call la_zlacgv(q - i + 1,x21(i,i),ldx21)
              call la_zlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              c = real(x21(i,i),KIND=dp)
              x21(i,i) = cone
              call la_zlarf('R',p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_zlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              call la_zlacgv(q - i + 1,x21(i,i),ldx21)
              if (i < m - q) then
                 s = sqrt(la_dznrm2(p - i,x11(i + 1,i),1)**2 + la_dznrm2(m - p - i,x21(i + 1, &
                           i),1)**2)
                 phi(i) = atan2(s,c)
              end if
           end do
           ! reduce the bottom-right portion of x11 to [ i 0 ]
           do i = m - q + 1,p
              call la_zlacgv(q - i + 1,x11(i,i),ldx11)
              call la_zlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              x11(i,i) = cone
              call la_zlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_zlarf('R',q - p,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(m - q + 1,i),ldx21, &
                        work(ilarf))
              call la_zlacgv(q - i + 1,x11(i,i),ldx11)
           end do
           ! reduce the bottom-right portion of x21 to [ 0 i ]
           do i = p + 1,q
              call la_zlacgv(q - i + 1,x21(m - q + i - p,i),ldx21)
              call la_zlarfgp(q - i + 1,x21(m - q + i - p,i),x21(m - q + i - p,i + 1),ldx21,tauq1(i))

              x21(m - q + i - p,i) = cone
              call la_zlarf('R',q - i,q - i + 1,x21(m - q + i - p,i),ldx21,tauq1(i),x21(m - q + i - p + 1,i) &
                        ,ldx21,work(ilarf))
              call la_zlacgv(q - i + 1,x21(m - q + i - p,i),ldx21)
           end do
           return
     end subroutine la_zunbdb4
     !> WUNBDB4: simultaneously bidiagonalizes the blocks of a tall and skinny
     !> matrix X with orthonomal columns:
     !> [ B11 ]
     !> [ X11 ]   [ P1 |    ] [  0  ]
     !> [-----] = [---------] [-----] Q1**T .
     !> [ X21 ]   [    | P2 ] [ B21 ]
     !> [  0  ]
     !> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
     !> M-P, or Q. Routines WUNBDB1, WUNBDB2, and WUNBDB3 handle cases in
     !> which M-Q is not the minimum dimension.
     !> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
     !> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
     !> Householder vectors.
     !> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
     !> implicitly by angles THETA, PHI.

     subroutine la_wunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,phi,taup1,taup2,tauq1, &
               phantom,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lwork,m,p,q,ldx11,ldx21
           ! Array Arguments
           real(qp),intent(out) :: phi(*),theta(*)
           complex(qp),intent(out) :: phantom(*),taup1(*),taup2(*),tauq1(*),work(*)
           complex(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
        ! ====================================================================

           ! Local Scalars
           real(qp) :: c,s
           integer(ilp) :: childinfo,i,ilarf,iorbdb5,j,llarf,lorbdb5,lworkmin, &
                     lworkopt
           logical(lk) :: lquery
           ! Intrinsic Function
           intrinsic :: atan2,cos,max,sin,sqrt
           ! Executable Statements
           ! test input arguments
           info = 0
           lquery = lwork == -1
           if (m < 0) then
              info = -1
           else if (p < m - q .or. m - p < m - q) then
              info = -2
           else if (q < m - q .or. q > m) then
              info = -3
           else if (ldx11 < max(1,p)) then
              info = -5
           else if (ldx21 < max(1,m - p)) then
              info = -7
           end if
           ! compute workspace
           if (info == 0) then
              ilarf = 2
              llarf = max(q - 1,p - 1,m - p - 1)
              iorbdb5 = 2
              lorbdb5 = q
              lworkopt = ilarf + llarf - 1
              lworkopt = max(lworkopt,iorbdb5 + lorbdb5 - 1)
              lworkmin = lworkopt
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WUNBDB4',-info)
              return
           else if (lquery) then
              return
           end if
           ! reduce columns 1, ..., m-q of x11 and x21
           do i = 1,m - q
              if (i == 1) then
                 do j = 1,m
                    phantom(j) = czero
                 end do
                 call la_wunbdb5(p,m - p,q,phantom(1),1,phantom(p + 1),1,x11,ldx11,x21, &
                           ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_wscal(p,cnegone,phantom(1),1)
                 call la_wlarfgp(p,phantom(1),phantom(2),1,taup1(1))
                 call la_wlarfgp(m - p,phantom(p + 1),phantom(p + 2),1,taup2(1))
                 theta(i) = atan2(real(phantom(1),KIND=qp),real(phantom(p + 1),KIND=qp))

                 c = cos(theta(i))
                 s = sin(theta(i))
                 phantom(1) = cone
                 phantom(p + 1) = cone
                 call la_wlarf('L',p,q,phantom(1),1,conjg(taup1(1)),x11,ldx11,work( &
                           ilarf))
                 call la_wlarf('L',m - p,q,phantom(p + 1),1,conjg(taup2(1)),x21,ldx21, &
                           work(ilarf))
              else
                 call la_wunbdb5(p - i + 1,m - p - i + 1,q - i + 1,x11(i,i - 1),1,x21(i,i - 1),1,x11(i,i) &
                           ,ldx11,x21(i,i),ldx21,work(iorbdb5),lorbdb5,childinfo)
                 call la_wscal(p - i + 1,cnegone,x11(i,i - 1),1)
                 call la_wlarfgp(p - i + 1,x11(i,i - 1),x11(i + 1,i - 1),1,taup1(i))
                 call la_wlarfgp(m - p - i + 1,x21(i,i - 1),x21(i + 1,i - 1),1,taup2(i))
                 theta(i) = atan2(real(x11(i,i - 1),KIND=qp),real(x21(i,i - 1),KIND=qp))
                 c = cos(theta(i))
                 s = sin(theta(i))
                 x11(i,i - 1) = cone
                 x21(i,i - 1) = cone
                 call la_wlarf('L',p - i + 1,q - i + 1,x11(i,i - 1),1,conjg(taup1(i)),x11(i,i), &
                           ldx11,work(ilarf))
                 call la_wlarf('L',m - p - i + 1,q - i + 1,x21(i,i - 1),1,conjg(taup2(i)),x21(i,i), &
                           ldx21,work(ilarf))
              end if
              call la_wqrot(q - i + 1,x11(i,i),ldx11,x21(i,i),ldx21,s,-c)
              call la_wlacgv(q - i + 1,x21(i,i),ldx21)
              call la_wlarfgp(q - i + 1,x21(i,i),x21(i,i + 1),ldx21,tauq1(i))
              c = real(x21(i,i),KIND=qp)
              x21(i,i) = cone
              call la_wlarf('R',p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_wlarf('R',m - p - i,q - i + 1,x21(i,i),ldx21,tauq1(i),x21(i + 1,i),ldx21, &
                        work(ilarf))
              call la_wlacgv(q - i + 1,x21(i,i),ldx21)
              if (i < m - q) then
                 s = sqrt(la_qwnrm2(p - i,x11(i + 1,i),1)**2 + la_qwnrm2(m - p - i,x21(i + 1, &
                           i),1)**2)
                 phi(i) = atan2(s,c)
              end if
           end do
           ! reduce the bottom-right portion of x11 to [ i 0 ]
           do i = m - q + 1,p
              call la_wlacgv(q - i + 1,x11(i,i),ldx11)
              call la_wlarfgp(q - i + 1,x11(i,i),x11(i,i + 1),ldx11,tauq1(i))
              x11(i,i) = cone
              call la_wlarf('R',p - i,q - i + 1,x11(i,i),ldx11,tauq1(i),x11(i + 1,i),ldx11, &
                        work(ilarf))
              call la_wlarf('R',q - p,q - i + 1,x11(i,i),ldx11,tauq1(i),x21(m - q + 1,i),ldx21, &
                        work(ilarf))
              call la_wlacgv(q - i + 1,x11(i,i),ldx11)
           end do
           ! reduce the bottom-right portion of x21 to [ 0 i ]
           do i = p + 1,q
              call la_wlacgv(q - i + 1,x21(m - q + i - p,i),ldx21)
              call la_wlarfgp(q - i + 1,x21(m - q + i - p,i),x21(m - q + i - p,i + 1),ldx21,tauq1(i))

              x21(m - q + i - p,i) = cone
              call la_wlarf('R',q - i,q - i + 1,x21(m - q + i - p,i),ldx21,tauq1(i),x21(m - q + i - p + 1,i) &
                        ,ldx21,work(ilarf))
              call la_wlacgv(q - i + 1,x21(m - q + i - p,i),ldx21)
           end do
           return
     end subroutine la_wunbdb4

     !> CUNCSD2BY1: computes the CS decomposition of an M-by-Q matrix X with
     !> orthonormal columns that has been partitioned into a 2-by-1 block
     !> structure:
     !> [  I1 0  0 ]
     !> [  0  C  0 ]
     !> [ X11 ]   [ U1 |    ] [  0  0  0 ]
     !> X = [-----] = [---------] [----------] V1**T .
     !> [ X21 ]   [    | U2 ] [  0  0  0 ]
     !> [  0  S  0 ]
     !> [  0  0  I2]
     !> X11 is P-by-Q. The unitary matrices U1, U2, and V1 are P-by-P,
     !> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     !> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     !> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     !> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).

     subroutine la_cuncsd2by1(jobu1,jobu2,jobv1t,m,p,q,x11,ldx11,x21,ldx21,theta, &
               u1,ldu1,u2,ldu2,v1t,ldv1t,work,lwork,rwork,lrwork,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,lwork,ldx11,ldx21,m,p,q
           integer(ilp),intent(in) :: lrwork
           integer(ilp) :: lrworkmin,lrworkopt
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           real(sp),intent(out) :: theta(*)
           complex(sp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),work(*)
           complex(sp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
           integer(ilp),intent(out) :: iwork(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,j,lbbcsd,lorbdb, &
           lorglq,lorglqmin,lorglqopt,lorgqr,lorgqrmin,lorgqropt,lworkmin,lworkopt, &
                     r
           logical(lk) :: lquery,wantu1,wantu2,wantv1t
           ! Local Arrays
           real(sp) :: dum(1)
           complex(sp) :: cdum(1,1)
           ! Intrinsic Function
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           lquery = (lwork == -1) .or. (lrwork == -1)
           if (m < 0) then
              info = -4
           else if (p < 0 .or. p > m) then
              info = -5
           else if (q < 0 .or. q > m) then
              info = -6
           else if (ldx11 < max(1,p)) then
              info = -8
           else if (ldx21 < max(1,m - p)) then
              info = -10
           else if (wantu1 .and. ldu1 < max(1,p)) then
              info = -13
           else if (wantu2 .and. ldu2 < max(1,m - p)) then
              info = -15
           else if (wantv1t .and. ldv1t < max(1,q)) then
              info = -17
           end if
           r = min(p,m - p,q,m - q)
           ! compute workspace
             ! work layout:
           ! |-----------------------------------------|
           ! | lworkopt (1)                            |
           ! |-----------------------------------------|
           ! | taup1 (max(1,p))                        |
           ! | taup2 (max(1,m-p))                      |
           ! | tauq1 (max(1,q))                        |
           ! |-----------------------------------------|
           ! | la_cunbdb work | la_cungqr work | la_cunglq work |
           ! |             |             |             |
           ! |             |             |             |
           ! |             |             |             |
           ! |             |             |             |
           ! |-----------------------------------------|
             ! rwork layout:
           ! |------------------|
           ! | lrworkopt (1)    |
           ! |------------------|
           ! | phi (max(1,r-1)) |
           ! |------------------|
           ! | b11d (r)         |
           ! | b11e (r-1)       |
           ! | b12d (r)         |
           ! | b12e (r-1)       |
           ! | b21d (r)         |
           ! | b21e (r-1)       |
           ! | b22d (r)         |
           ! | b22e (r-1)       |
           ! | la_cbbcsd rwork     |
           ! |------------------|
           if (info == 0) then
              iphi = 2
              ib11d = iphi + max(1,r - 1)
              ib11e = ib11d + max(1,r)
              ib12d = ib11e + max(1,r - 1)
              ib12e = ib12d + max(1,r)
              ib21d = ib12e + max(1,r - 1)
              ib21e = ib21d + max(1,r)
              ib22d = ib21e + max(1,r - 1)
              ib22e = ib22d + max(1,r)
              ibbcsd = ib22e + max(1,r - 1)
              itaup1 = 2
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              iorbdb = itauq1 + max(1,q)
              iorgqr = itauq1 + max(1,q)
              iorglq = itauq1 + max(1,q)
              lorgqrmin = 1
              lorgqropt = 1
              lorglqmin = 1
              lorglqopt = 1
              if (r == q) then
                 call la_cunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work,-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_cungqr(p,p,q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_cungqr(m - p,m - p,q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_cunglq(q - 1,q - 1,q - 1,v1t,ldv1t,cdum,work(1),-1,childinfo)

                    lorglqmin = max(lorglqmin,q - 1)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_cbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,dum(1),u1, &
                 ldu1,u2,ldu2,v1t,ldv1t,cdum,1,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else if (r == p) then
                 call la_cunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_cungqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,cdum,work(1),-1,childinfo &
                              )
                    lorgqrmin = max(lorgqrmin,p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_cungqr(m - p,m - p,q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_cunglq(q,q,r,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_cbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,dum,v1t, &
                 ldv1t,cdum,1,u1,ldu1,u2,ldu2,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else if (r == m - p) then
                 call la_cunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_cungqr(p,p,q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_cungqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,cdum,work(1),-1, &
                              childinfo)
                    lorgqrmin = max(lorgqrmin,m - p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_cunglq(q,q,r,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_cbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,dum,cdum, &
                  1,v1t,ldv1t,u2,ldu2,u1,ldu1,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                            1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else
                 call la_cunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,cdum,work(1),-1,childinfo)
                 lorbdb = m + int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_cungqr(p,p,m - q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_cungqr(m - p,m - p,m - q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_cunglq(q,q,q,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_cbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,dum,u2, &
                 ldu2,u1,ldu1,cdum,1,v1t,ldv1t,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              end if
              lrworkmin = ibbcsd + lbbcsd - 1
              lrworkopt = lrworkmin
              rwork(1) = lrworkopt
              lworkmin = max(iorbdb + lorbdb - 1,iorgqr + lorgqrmin - 1,iorglq + lorglqmin - 1)
              lworkopt = max(iorbdb + lorbdb - 1,iorgqr + lorgqropt - 1,iorglq + lorglqopt - 1)
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -19
              end if
              if (lrwork < lrworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CUNCSD2BY1',-info)
              return
           else if (lquery) then
              return
           end if
           lorgqr = lwork - iorgqr + 1
           lorglq = lwork - iorglq + 1
           ! handle four cases separately: r = q, r = p, r = m-p, and r = m-q,
           ! in which r = min(p,m-p,q,m-q)
           if (r == q) then
              ! case 1: r = q
              ! simultaneously bidiagonalize x11 and x21
              call la_cunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_clacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_cungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_clacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_cungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_clacpy('U',q - 1,q - 1,x21(1,2),ldx21,v1t(2,2),ldv1t)
                 call la_cunglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglq,childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_cbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,rwork(iphi),u1, &
              ldu1,u2,ldu2,v1t,ldv1t,cdum,1,rwork(ib11d),rwork(ib11e),rwork(ib12d),rwork( &
              ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd),lrwork - &
                        ibbcsd + 1,childinfo)
              ! permute rows and columns to place czero submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_clapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == p) then
              ! case 2: r = p
              ! simultaneously bidiagonalize x11 and x21
              call la_cunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 u1(1,1) = cone
                 do j = 2,p
                    u1(1,j) = czero
                    u1(j,1) = czero
                 end do
                 call la_clacpy('L',p - 1,p - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_cungqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,work(itaup1),work(iorgqr), &
                           lorgqr,childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_clacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_cungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_clacpy('U',p,q,x11,ldx11,v1t,ldv1t)
                 call la_cunglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_cbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,rwork(iphi),v1t, &
               ldv1t,cdum,1,u1,ldu1,u2,ldu2,rwork(ib11d),rwork(ib11e),rwork(ib12d),rwork( &
               ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                         lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_clapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == m - p) then
              ! case 3: r = m-p
              ! simultaneously bidiagonalize x11 and x21
              call la_cunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_clacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_cungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 u2(1,1) = cone
                 do j = 2,m - p
                    u2(1,j) = czero
                    u2(j,1) = czero
                 end do
                 call la_clacpy('L',m - p - 1,m - p - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_cungqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,work(itaup2),work(iorgqr) &
                           ,lorgqr,childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_clacpy('U',m - p,q,x21,ldx21,v1t,ldv1t)
                 call la_cunglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_cbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,rwork(iphi), &
              cdum,1,v1t,ldv1t,u2,ldu2,u1,ldu1,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
              rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                        lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > r) then
                 do i = 1,r
                    iwork(i) = q - r + i
                 end do
                 do i = r + 1,q
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_clapmt(.false.,p,q,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_clapmr(.false.,q,q,v1t,ldv1t,iwork)
                 end if
              end if
           else
              ! case 4: r = m-q
              ! simultaneously bidiagonalize x11 and x21
              call la_cunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
              itaup1),work(itaup2),work(itauq1),work(iorbdb),work(iorbdb + m),lorbdb - m, &
                        childinfo)
              ! accumulate householder reflectors
              if (wantu2 .and. m - p > 0) then
                 call la_ccopy(m - p,work(iorbdb + p),1,u2,1)
              end if
              if (wantu1 .and. p > 0) then
                 call la_ccopy(p,work(iorbdb),1,u1,1)
                 do j = 2,p
                    u1(1,j) = czero
                 end do
                 call la_clacpy('L',p - 1,m - q - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_cungqr(p,p,m - q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 do j = 2,m - p
                    u2(1,j) = czero
                 end do
                 call la_clacpy('L',m - p - 1,m - q - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_cungqr(m - p,m - p,m - q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_clacpy('U',m - q,q,x21,ldx21,v1t,ldv1t)
                 call la_clacpy('U',p - (m - q),q - (m - q),x11(m - q + 1,m - q + 1),ldx11,v1t(m - q + 1,m - q + &
                           1),ldv1t)
                 call la_clacpy('U',-p + q,q - p,x21(m - q + 1,p + 1),ldx21,v1t(p + 1,p + 1),ldv1t)

                 call la_cunglq(q,q,q,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_cbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,rwork(iphi), &
              u2,ldu2,u1,ldu1,cdum,1,v1t,ldv1t,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
              rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                        lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (p > r) then
                 do i = 1,r
                    iwork(i) = p - r + i
                 end do
                 do i = r + 1,p
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_clapmt(.false.,p,p,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_clapmr(.false.,p,q,v1t,ldv1t,iwork)
                 end if
              end if
           end if
           return
     end subroutine la_cuncsd2by1
     !> ZUNCSD2BY1: computes the CS decomposition of an M-by-Q matrix X with
     !> orthonormal columns that has been partitioned into a 2-by-1 block
     !> structure:
     !> [  I1 0  0 ]
     !> [  0  C  0 ]
     !> [ X11 ]   [ U1 |    ] [  0  0  0 ]
     !> X = [-----] = [---------] [----------] V1**T .
     !> [ X21 ]   [    | U2 ] [  0  0  0 ]
     !> [  0  S  0 ]
     !> [  0  0  I2]
     !> X11 is P-by-Q. The unitary matrices U1, U2, and V1 are P-by-P,
     !> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     !> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     !> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     !> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).

     subroutine la_zuncsd2by1(jobu1,jobu2,jobv1t,m,p,q,x11,ldx11,x21,ldx21,theta, &
               u1,ldu1,u2,ldu2,v1t,ldv1t,work,lwork,rwork,lrwork,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,lwork,ldx11,ldx21,m,p,q
           integer(ilp),intent(in) :: lrwork
           integer(ilp) :: lrworkmin,lrworkopt
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           real(dp),intent(out) :: theta(*)
           complex(dp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),work(*)
           complex(dp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
           integer(ilp),intent(out) :: iwork(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,j,lbbcsd,lorbdb, &
           lorglq,lorglqmin,lorglqopt,lorgqr,lorgqrmin,lorgqropt,lworkmin,lworkopt, &
                     r
           logical(lk) :: lquery,wantu1,wantu2,wantv1t
           ! Local Arrays
           real(dp) :: dum(1)
           complex(dp) :: cdum(1,1)
           ! Intrinsic Function
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           lquery = (lwork == -1) .or. (lrwork == -1)
           if (m < 0) then
              info = -4
           else if (p < 0 .or. p > m) then
              info = -5
           else if (q < 0 .or. q > m) then
              info = -6
           else if (ldx11 < max(1,p)) then
              info = -8
           else if (ldx21 < max(1,m - p)) then
              info = -10
           else if (wantu1 .and. ldu1 < max(1,p)) then
              info = -13
           else if (wantu2 .and. ldu2 < max(1,m - p)) then
              info = -15
           else if (wantv1t .and. ldv1t < max(1,q)) then
              info = -17
           end if
           r = min(p,m - p,q,m - q)
           ! compute workspace
             ! work layout:
           ! |-----------------------------------------|
           ! | lworkopt (1)                            |
           ! |-----------------------------------------|
           ! | taup1 (max(1,p))                        |
           ! | taup2 (max(1,m-p))                      |
           ! | tauq1 (max(1,q))                        |
           ! |-----------------------------------------|
           ! | la_zunbdb work | la_zungqr work | la_zunglq work |
           ! |             |             |             |
           ! |             |             |             |
           ! |             |             |             |
           ! |             |             |             |
           ! |-----------------------------------------|
             ! rwork layout:
           ! |------------------|
           ! | lrworkopt (1)    |
           ! |------------------|
           ! | phi (max(1,r-1)) |
           ! |------------------|
           ! | b11d (r)         |
           ! | b11e (r-1)       |
           ! | b12d (r)         |
           ! | b12e (r-1)       |
           ! | b21d (r)         |
           ! | b21e (r-1)       |
           ! | b22d (r)         |
           ! | b22e (r-1)       |
           ! | la_zbbcsd rwork     |
           ! |------------------|
           if (info == 0) then
              iphi = 2
              ib11d = iphi + max(1,r - 1)
              ib11e = ib11d + max(1,r)
              ib12d = ib11e + max(1,r - 1)
              ib12e = ib12d + max(1,r)
              ib21d = ib12e + max(1,r - 1)
              ib21e = ib21d + max(1,r)
              ib22d = ib21e + max(1,r - 1)
              ib22e = ib22d + max(1,r)
              ibbcsd = ib22e + max(1,r - 1)
              itaup1 = 2
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              iorbdb = itauq1 + max(1,q)
              iorgqr = itauq1 + max(1,q)
              iorglq = itauq1 + max(1,q)
              lorgqrmin = 1
              lorgqropt = 1
              lorglqmin = 1
              lorglqopt = 1
              if (r == q) then
                 call la_zunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work,-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_zungqr(p,p,q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_zungqr(m - p,m - p,q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_zunglq(q - 1,q - 1,q - 1,v1t,ldv1t,cdum,work(1),-1,childinfo)

                    lorglqmin = max(lorglqmin,q - 1)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_zbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,dum,u1,ldu1, &
                  u2,ldu2,v1t,ldv1t,cdum,1,dum,dum,dum,dum,dum,dum,dum,dum,rwork(1),- &
                            1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else if (r == p) then
                 call la_zunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_zungqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,cdum,work(1),-1,childinfo &
                              )
                    lorgqrmin = max(lorgqrmin,p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_zungqr(m - p,m - p,q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_zunglq(q,q,r,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_zbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,dum,v1t, &
                 ldv1t,cdum,1,u1,ldu1,u2,ldu2,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else if (r == m - p) then
                 call la_zunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_zungqr(p,p,q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_zungqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,cdum,work(1),-1, &
                              childinfo)
                    lorgqrmin = max(lorgqrmin,m - p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_zunglq(q,q,r,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_zbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,dum,cdum, &
                  1,v1t,ldv1t,u2,ldu2,u1,ldu1,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                            1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else
                 call la_zunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,cdum,work(1),-1,childinfo)
                 lorbdb = m + int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_zungqr(p,p,m - q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_zungqr(m - p,m - p,m - q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_zunglq(q,q,q,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_zbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,dum,u2, &
                 ldu2,u1,ldu1,cdum,1,v1t,ldv1t,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              end if
              lrworkmin = ibbcsd + lbbcsd - 1
              lrworkopt = lrworkmin
              rwork(1) = lrworkopt
              lworkmin = max(iorbdb + lorbdb - 1,iorgqr + lorgqrmin - 1,iorglq + lorglqmin - 1)
              lworkopt = max(iorbdb + lorbdb - 1,iorgqr + lorgqropt - 1,iorglq + lorglqopt - 1)
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -19
              end if
              if (lrwork < lrworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZUNCSD2BY1',-info)
              return
           else if (lquery) then
              return
           end if
           lorgqr = lwork - iorgqr + 1
           lorglq = lwork - iorglq + 1
           ! handle four cases separately: r = q, r = p, r = m-p, and r = m-q,
           ! in which r = min(p,m-p,q,m-q)
           if (r == q) then
              ! case 1: r = q
              ! simultaneously bidiagonalize x11 and x21
              call la_zunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_zlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_zungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_zlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_zungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_zlacpy('U',q - 1,q - 1,x21(1,2),ldx21,v1t(2,2),ldv1t)
                 call la_zunglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglq,childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_zbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,rwork(iphi),u1, &
              ldu1,u2,ldu2,v1t,ldv1t,cdum,1,rwork(ib11d),rwork(ib11e),rwork(ib12d),rwork( &
              ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd),lrwork - &
                        ibbcsd + 1,childinfo)
              ! permute rows and columns to place czero submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_zlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == p) then
              ! case 2: r = p
              ! simultaneously bidiagonalize x11 and x21
              call la_zunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 u1(1,1) = cone
                 do j = 2,p
                    u1(1,j) = czero
                    u1(j,1) = czero
                 end do
                 call la_zlacpy('L',p - 1,p - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_zungqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,work(itaup1),work(iorgqr), &
                           lorgqr,childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_zlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_zungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_zlacpy('U',p,q,x11,ldx11,v1t,ldv1t)
                 call la_zunglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_zbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,rwork(iphi),v1t, &
               ldv1t,cdum,1,u1,ldu1,u2,ldu2,rwork(ib11d),rwork(ib11e),rwork(ib12d),rwork( &
               ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                         lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_zlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == m - p) then
              ! case 3: r = m-p
              ! simultaneously bidiagonalize x11 and x21
              call la_zunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_zlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_zungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 u2(1,1) = cone
                 do j = 2,m - p
                    u2(1,j) = czero
                    u2(j,1) = czero
                 end do
                 call la_zlacpy('L',m - p - 1,m - p - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_zungqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,work(itaup2),work(iorgqr) &
                           ,lorgqr,childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_zlacpy('U',m - p,q,x21,ldx21,v1t,ldv1t)
                 call la_zunglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_zbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,rwork(iphi), &
              cdum,1,v1t,ldv1t,u2,ldu2,u1,ldu1,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
              rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                        lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > r) then
                 do i = 1,r
                    iwork(i) = q - r + i
                 end do
                 do i = r + 1,q
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_zlapmt(.false.,p,q,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_zlapmr(.false.,q,q,v1t,ldv1t,iwork)
                 end if
              end if
           else
              ! case 4: r = m-q
              ! simultaneously bidiagonalize x11 and x21
              call la_zunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
              itaup1),work(itaup2),work(itauq1),work(iorbdb),work(iorbdb + m),lorbdb - m, &
                        childinfo)
              ! accumulate householder reflectors
              if (wantu2 .and. m - p > 0) then
                 call la_zcopy(m - p,work(iorbdb + p),1,u2,1)
              end if
              if (wantu1 .and. p > 0) then
                 call la_zcopy(p,work(iorbdb),1,u1,1)
                 do j = 2,p
                    u1(1,j) = czero
                 end do
                 call la_zlacpy('L',p - 1,m - q - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_zungqr(p,p,m - q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 do j = 2,m - p
                    u2(1,j) = czero
                 end do
                 call la_zlacpy('L',m - p - 1,m - q - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_zungqr(m - p,m - p,m - q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_zlacpy('U',m - q,q,x21,ldx21,v1t,ldv1t)
                 call la_zlacpy('U',p - (m - q),q - (m - q),x11(m - q + 1,m - q + 1),ldx11,v1t(m - q + 1,m - q + &
                           1),ldv1t)
                 call la_zlacpy('U',-p + q,q - p,x21(m - q + 1,p + 1),ldx21,v1t(p + 1,p + 1),ldv1t)

                 call la_zunglq(q,q,q,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_zbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,rwork(iphi), &
              u2,ldu2,u1,ldu1,cdum,1,v1t,ldv1t,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
              rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                        lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (p > r) then
                 do i = 1,r
                    iwork(i) = p - r + i
                 end do
                 do i = r + 1,p
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_zlapmt(.false.,p,p,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_zlapmr(.false.,p,q,v1t,ldv1t,iwork)
                 end if
              end if
           end if
           return
     end subroutine la_zuncsd2by1
     !> WUNCSD2BY1: computes the CS decomposition of an M-by-Q matrix X with
     !> orthonormal columns that has been partitioned into a 2-by-1 block
     !> structure:
     !> [  I1 0  0 ]
     !> [  0  C  0 ]
     !> [ X11 ]   [ U1 |    ] [  0  0  0 ]
     !> X = [-----] = [---------] [----------] V1**T .
     !> [ X21 ]   [    | U2 ] [  0  0  0 ]
     !> [  0  S  0 ]
     !> [  0  0  I2]
     !> X11 is P-by-Q. The unitary matrices U1, U2, and V1 are P-by-P,
     !> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
     !> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
     !> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
     !> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).

     subroutine la_wuncsd2by1(jobu1,jobu2,jobv1t,m,p,q,x11,ldx11,x21,ldx21,theta, &
               u1,ldu1,u2,ldu2,v1t,ldv1t,work,lwork,rwork,lrwork,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobu1,jobu2,jobv1t
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu1,ldu2,ldv1t,lwork,ldx11,ldx21,m,p,q
           integer(ilp),intent(in) :: lrwork
           integer(ilp) :: lrworkmin,lrworkopt
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           real(qp),intent(out) :: theta(*)
           complex(qp),intent(out) :: u1(ldu1,*),u2(ldu2,*),v1t(ldv1t,*),work(*)
           complex(qp),intent(inout) :: x11(ldx11,*),x21(ldx21,*)
           integer(ilp),intent(out) :: iwork(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: childinfo,i,ib11d,ib11e,ib12d,ib12e,ib21d,ib21e,ib22d,ib22e, &
           ibbcsd,iorbdb,iorglq,iorgqr,iphi,itaup1,itaup2,itauq1,j,lbbcsd,lorbdb, &
           lorglq,lorglqmin,lorglqopt,lorgqr,lorgqrmin,lorgqropt,lworkmin,lworkopt, &
                     r
           logical(lk) :: lquery,wantu1,wantu2,wantv1t
           ! Local Arrays
           real(qp) :: dum(1)
           complex(qp) :: cdum(1,1)
           ! Intrinsic Function
           intrinsic :: int,max,min
           ! Executable Statements
           ! test input arguments
           info = 0
           wantu1 = la_lsame(jobu1,'Y')
           wantu2 = la_lsame(jobu2,'Y')
           wantv1t = la_lsame(jobv1t,'Y')
           lquery = (lwork == -1) .or. (lrwork == -1)
           if (m < 0) then
              info = -4
           else if (p < 0 .or. p > m) then
              info = -5
           else if (q < 0 .or. q > m) then
              info = -6
           else if (ldx11 < max(1,p)) then
              info = -8
           else if (ldx21 < max(1,m - p)) then
              info = -10
           else if (wantu1 .and. ldu1 < max(1,p)) then
              info = -13
           else if (wantu2 .and. ldu2 < max(1,m - p)) then
              info = -15
           else if (wantv1t .and. ldv1t < max(1,q)) then
              info = -17
           end if
           r = min(p,m - p,q,m - q)
           ! compute workspace
             ! work layout:
           ! |-----------------------------------------|
           ! | lworkopt (1)                            |
           ! |-----------------------------------------|
           ! | taup1 (max(1,p))                        |
           ! | taup2 (max(1,m-p))                      |
           ! | tauq1 (max(1,q))                        |
           ! |-----------------------------------------|
           ! | la_wunbdb work | la_wungqr work | la_wunglq work |
           ! |             |             |             |
           ! |             |             |             |
           ! |             |             |             |
           ! |             |             |             |
           ! |-----------------------------------------|
             ! rwork layout:
           ! |------------------|
           ! | lrworkopt (1)    |
           ! |------------------|
           ! | phi (max(1,r-1)) |
           ! |------------------|
           ! | b11d (r)         |
           ! | b11e (r-1)       |
           ! | b12d (r)         |
           ! | b12e (r-1)       |
           ! | b21d (r)         |
           ! | b21e (r-1)       |
           ! | b22d (r)         |
           ! | b22e (r-1)       |
           ! | la_wbbcsd rwork     |
           ! |------------------|
           if (info == 0) then
              iphi = 2
              ib11d = iphi + max(1,r - 1)
              ib11e = ib11d + max(1,r)
              ib12d = ib11e + max(1,r - 1)
              ib12e = ib12d + max(1,r)
              ib21d = ib12e + max(1,r - 1)
              ib21e = ib21d + max(1,r)
              ib22d = ib21e + max(1,r - 1)
              ib22e = ib22d + max(1,r)
              ibbcsd = ib22e + max(1,r - 1)
              itaup1 = 2
              itaup2 = itaup1 + max(1,p)
              itauq1 = itaup2 + max(1,m - p)
              iorbdb = itauq1 + max(1,q)
              iorgqr = itauq1 + max(1,q)
              iorglq = itauq1 + max(1,q)
              lorgqrmin = 1
              lorgqropt = 1
              lorglqmin = 1
              lorglqopt = 1
              if (r == q) then
                 call la_wunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work,-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_wungqr(p,p,q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_wungqr(m - p,m - p,q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_wunglq(q - 1,q - 1,q - 1,v1t,ldv1t,cdum,work(1),-1,childinfo)

                    lorglqmin = max(lorglqmin,q - 1)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_wbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,dum,u1,ldu1, &
                  u2,ldu2,v1t,ldv1t,cdum,1,dum,dum,dum,dum,dum,dum,dum,dum,rwork(1),- &
                            1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else if (r == p) then
                 call la_wunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_wungqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,cdum,work(1),-1,childinfo &
                              )
                    lorgqrmin = max(lorgqrmin,p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_wungqr(m - p,m - p,q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_wunglq(q,q,r,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_wbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,dum,v1t, &
                 ldv1t,cdum,1,u1,ldu1,u2,ldu2,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else if (r == m - p) then
                 call la_wunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,work(1),-1,childinfo)
                 lorbdb = int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_wungqr(p,p,q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_wungqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,cdum,work(1),-1, &
                              childinfo)
                    lorgqrmin = max(lorgqrmin,m - p - 1)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_wunglq(q,q,r,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_wbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,dum,cdum, &
                  1,v1t,ldv1t,u2,ldu2,u1,ldu1,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                            1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              else
                 call la_wunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,dum,cdum,cdum, &
                           cdum,cdum,work(1),-1,childinfo)
                 lorbdb = m + int(work(1),KIND=ilp)
                 if (wantu1 .and. p > 0) then
                    call la_wungqr(p,p,m - q,u1,ldu1,cdum,work(1),-1,childinfo)
                    lorgqrmin = max(lorgqrmin,p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantu2 .and. m - p > 0) then
                    call la_wungqr(m - p,m - p,m - q,u2,ldu2,cdum,work(1),-1,childinfo)

                    lorgqrmin = max(lorgqrmin,m - p)
                    lorgqropt = max(lorgqropt,int(work(1),KIND=ilp))
                 end if
                 if (wantv1t .and. q > 0) then
                    call la_wunglq(q,q,q,v1t,ldv1t,cdum,work(1),-1,childinfo)
                    lorglqmin = max(lorglqmin,q)
                    lorglqopt = max(lorglqopt,int(work(1),KIND=ilp))
                 end if
                 call la_wbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,dum,u2, &
                 ldu2,u1,ldu1,cdum,1,v1t,ldv1t,dum,dum,dum,dum,dum,dum,dum,dum,rwork( &
                           1),-1,childinfo)
                 lbbcsd = int(rwork(1),KIND=ilp)
              end if
              lrworkmin = ibbcsd + lbbcsd - 1
              lrworkopt = lrworkmin
              rwork(1) = lrworkopt
              lworkmin = max(iorbdb + lorbdb - 1,iorgqr + lorgqrmin - 1,iorglq + lorglqmin - 1)
              lworkopt = max(iorbdb + lorbdb - 1,iorgqr + lorgqropt - 1,iorglq + lorglqopt - 1)
              work(1) = lworkopt
              if (lwork < lworkmin .and. .not. lquery) then
                 info = -19
              end if
              if (lrwork < lrworkmin .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WUNCSD2BY1',-info)
              return
           else if (lquery) then
              return
           end if
           lorgqr = lwork - iorgqr + 1
           lorglq = lwork - iorglq + 1
           ! handle four cases separately: r = q, r = p, r = m-p, and r = m-q,
           ! in which r = min(p,m-p,q,m-q)
           if (r == q) then
              ! case 1: r = q
              ! simultaneously bidiagonalize x11 and x21
              call la_wunbdb1(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_wlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_wungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_wlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_wungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 v1t(1,1) = cone
                 do j = 2,q
                    v1t(1,j) = czero
                    v1t(j,1) = czero
                 end do
                 call la_wlacpy('U',q - 1,q - 1,x21(1,2),ldx21,v1t(2,2),ldv1t)
                 call la_wunglq(q - 1,q - 1,q - 1,v1t(2,2),ldv1t,work(itauq1),work(iorglq), &
                           lorglq,childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_wbbcsd(jobu1,jobu2,jobv1t,'N','N',m,p,q,theta,rwork(iphi),u1, &
              ldu1,u2,ldu2,v1t,ldv1t,cdum,1,rwork(ib11d),rwork(ib11e),rwork(ib12d),rwork( &
              ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd),lrwork - &
                        ibbcsd + 1,childinfo)
              ! permute rows and columns to place czero submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_wlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == p) then
              ! case 2: r = p
              ! simultaneously bidiagonalize x11 and x21
              call la_wunbdb2(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 u1(1,1) = cone
                 do j = 2,p
                    u1(1,j) = czero
                    u1(j,1) = czero
                 end do
                 call la_wlacpy('L',p - 1,p - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_wungqr(p - 1,p - 1,p - 1,u1(2,2),ldu1,work(itaup1),work(iorgqr), &
                           lorgqr,childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 call la_wlacpy('L',m - p,q,x21,ldx21,u2,ldu2)
                 call la_wungqr(m - p,m - p,q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_wlacpy('U',p,q,x11,ldx11,v1t,ldv1t)
                 call la_wunglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_wbbcsd(jobv1t,'N',jobu1,jobu2,'T',m,q,p,theta,rwork(iphi),v1t, &
               ldv1t,cdum,1,u1,ldu1,u2,ldu2,rwork(ib11d),rwork(ib11e),rwork(ib12d),rwork( &
               ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                         lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > 0 .and. wantu2) then
                 do i = 1,q
                    iwork(i) = m - p - q + i
                 end do
                 do i = q + 1,m - p
                    iwork(i) = i - q
                 end do
                 call la_wlapmt(.false.,m - p,m - p,u2,ldu2,iwork)
              end if
           else if (r == m - p) then
              ! case 3: r = m-p
              ! simultaneously bidiagonalize x11 and x21
              call la_wunbdb3(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
                        itaup1),work(itaup2),work(itauq1),work(iorbdb),lorbdb,childinfo)
              ! accumulate householder reflectors
              if (wantu1 .and. p > 0) then
                 call la_wlacpy('L',p,q,x11,ldx11,u1,ldu1)
                 call la_wungqr(p,p,q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 u2(1,1) = cone
                 do j = 2,m - p
                    u2(1,j) = czero
                    u2(j,1) = czero
                 end do
                 call la_wlacpy('L',m - p - 1,m - p - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_wungqr(m - p - 1,m - p - 1,m - p - 1,u2(2,2),ldu2,work(itaup2),work(iorgqr) &
                           ,lorgqr,childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_wlacpy('U',m - p,q,x21,ldx21,v1t,ldv1t)
                 call la_wunglq(q,q,r,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_wbbcsd('N',jobv1t,jobu2,jobu1,'T',m,m - q,m - p,theta,rwork(iphi), &
              cdum,1,v1t,ldv1t,u2,ldu2,u1,ldu1,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
              rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                        lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (q > r) then
                 do i = 1,r
                    iwork(i) = q - r + i
                 end do
                 do i = r + 1,q
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_wlapmt(.false.,p,q,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_wlapmr(.false.,q,q,v1t,ldv1t,iwork)
                 end if
              end if
           else
              ! case 4: r = m-q
              ! simultaneously bidiagonalize x11 and x21
              call la_wunbdb4(m,p,q,x11,ldx11,x21,ldx21,theta,rwork(iphi),work( &
              itaup1),work(itaup2),work(itauq1),work(iorbdb),work(iorbdb + m),lorbdb - m, &
                        childinfo)
              ! accumulate householder reflectors
              if (wantu2 .and. m - p > 0) then
                 call la_wcopy(m - p,work(iorbdb + p),1,u2,1)
              end if
              if (wantu1 .and. p > 0) then
                 call la_wcopy(p,work(iorbdb),1,u1,1)
                 do j = 2,p
                    u1(1,j) = czero
                 end do
                 call la_wlacpy('L',p - 1,m - q - 1,x11(2,1),ldx11,u1(2,2),ldu1)
                 call la_wungqr(p,p,m - q,u1,ldu1,work(itaup1),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantu2 .and. m - p > 0) then
                 do j = 2,m - p
                    u2(1,j) = czero
                 end do
                 call la_wlacpy('L',m - p - 1,m - q - 1,x21(2,1),ldx21,u2(2,2),ldu2)
                 call la_wungqr(m - p,m - p,m - q,u2,ldu2,work(itaup2),work(iorgqr),lorgqr, &
                           childinfo)
              end if
              if (wantv1t .and. q > 0) then
                 call la_wlacpy('U',m - q,q,x21,ldx21,v1t,ldv1t)
                 call la_wlacpy('U',p - (m - q),q - (m - q),x11(m - q + 1,m - q + 1),ldx11,v1t(m - q + 1,m - q + &
                           1),ldv1t)
                 call la_wlacpy('U',-p + q,q - p,x21(m - q + 1,p + 1),ldx21,v1t(p + 1,p + 1),ldv1t)

                 call la_wunglq(q,q,q,v1t,ldv1t,work(itauq1),work(iorglq),lorglq, &
                           childinfo)
              end if
              ! simultaneously diagonalize x11 and x21.
              call la_wbbcsd(jobu2,jobu1,'N',jobv1t,'N',m,m - p,m - q,theta,rwork(iphi), &
              u2,ldu2,u1,ldu1,cdum,1,v1t,ldv1t,rwork(ib11d),rwork(ib11e),rwork(ib12d), &
              rwork(ib12e),rwork(ib21d),rwork(ib21e),rwork(ib22d),rwork(ib22e),rwork(ibbcsd), &
                        lbbcsd,childinfo)
              ! permute rows and columns to place identity submatrices in
              ! preferred positions
              if (p > r) then
                 do i = 1,r
                    iwork(i) = p - r + i
                 end do
                 do i = r + 1,p
                    iwork(i) = i - r
                 end do
                 if (wantu1) then
                    call la_wlapmt(.false.,p,p,u1,ldu1,iwork)
                 end if
                 if (wantv1t) then
                    call la_wlapmr(.false.,p,q,v1t,ldv1t,iwork)
                 end if
              end if
           end if
           return
     end subroutine la_wuncsd2by1

end module la_lapack_cosine_sine
